################################################################################
############### Grundliegend: Funktion, die aus einer langen Liste #############
################### eine Kontingenztabelle macht ###############################
################################################################################

make_contingency <- function(dataset,Project, Organization){
  
  network <- dataset
  names(network)<- c(as.character(Project),as.character(Organization))
  
  goalset<- list()
  
  for(i in 1:length(unique(network$Project))){
    dat<- network[network$Project==unique(network$Project)[i],]
    if(nrow(dat)==1){
      print(paste0("Projekte mit nur einem gemeinsamen Partner: ", unique(network$Project)[i]))
      dat2<- data.frame("V1" = dat$Org,
                        "V2" = dat$Org,
                        "Project" = dat$Project)             
      
    } else{
      dat2<- as.data.frame(combinations(n= length(unique(dat$Org)), r= 2 , v = dat$Org,repeats.allowed = F))
      dat2$Project<- NA
      dat2[,3]<- unique(network$Project)[i]
    }
    goalset[[i]]<- dat2
  }
  
  goalset_2<- do.call(rbind,goalset)
  
  return(goalset_2)
}


run_leiden_cluster_workflow <- function(
    graph,
    resolution_values = seq(0.2, 3, by = 0.05),
    n_runs = 100,
    min_clusters = 2,
    objective_function = "modularity",
    n_iterations = 10,
    seed = 123,
    use_weights = TRUE
) {
  
  if (!inherits(graph, "igraph")) {
    stop("`graph` must be an igraph object.")
  }
  
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package `igraph` is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required.")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("Package `tidyr` is required.")
  }
  if (!requireNamespace("purrr", quietly = TRUE)) {
    stop("Package `purrr` is required.")
  }
  
  set.seed(seed)
  
  edge_weights <- NULL
  
  if (use_weights) {
    if ("weight" %in% igraph::edge_attr_names(graph)) {
      edge_weights <- igraph::E(graph)$weight
    } else {
      warning("No edge attribute `weight` found. Running unweighted Leiden clustering.")
      edge_weights <- NULL
    }
  }
  
  calc_mean_nmi <- function(membership_list) {
    
    n <- length(membership_list)
    
    if (n < 2) return(NA_real_)
    
    combs <- utils::combn(seq_len(n), 2)
    
    nmi_values <- apply(combs, 2, function(idx) {
      igraph::compare(
        membership_list[[idx[1]]],
        membership_list[[idx[2]]],
        method = "nmi"
      )
    })
    
    mean(nmi_values, na.rm = TRUE)
  }
  
  leiden_runs <- tidyr::crossing(
    resolution = resolution_values,
    run = seq_len(n_runs)
  ) |>
    dplyr::mutate(
      cluster = purrr::map(
        resolution,
        \(res) igraph::cluster_leiden(
          graph,
          objective_function = objective_function,
          weights = edge_weights,
          resolution = res,
          n_iterations = n_iterations
        )
      ),
      membership = purrr::map(cluster, igraph::membership),
      n_clusters = purrr::map_int(
        membership,
        \(x) length(unique(x))
      ),
      modularity = purrr::map_dbl(
        membership,
        \(x) igraph::modularity(
          graph,
          membership = x,
          weights = edge_weights
        )
      )
    )
  
  leiden_summary <- leiden_runs |>
    dplyr::group_by(resolution, n_clusters) |>
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_modularity = mean(modularity, na.rm = TRUE),
      sd_modularity = stats::sd(modularity, na.rm = TRUE),
      max_modularity = max(modularity, na.rm = TRUE),
      .groups = "drop"
    )
  
  leiden_stability <- leiden_runs |>
    dplyr::group_by(resolution, n_clusters) |>
    dplyr::summarise(
      mean_nmi = calc_mean_nmi(membership),
      .groups = "drop"
    )
  
  leiden_eval <- leiden_summary |>
    dplyr::left_join(
      leiden_stability,
      by = c("resolution", "n_clusters")
    ) |>
    dplyr::arrange(
      dplyr::desc(mean_nmi),
      dplyr::desc(mean_modularity)
    )
  
  candidate_solutions <- leiden_eval |>
    dplyr::filter(n_clusters >= min_clusters) |>
    dplyr::arrange(
      dplyr::desc(mean_nmi),
      dplyr::desc(mean_modularity),
      n_clusters
    )
  
  if (nrow(candidate_solutions) == 0) {
    stop("No candidate solution found. Try lowering `min_clusters` or widening `resolution_values`.")
  }
  
  selected_resolution <- candidate_solutions$resolution[1]
  selected_n_clusters <- candidate_solutions$n_clusters[1]
  
  best_run <- leiden_runs |>
    dplyr::filter(
      resolution == selected_resolution,
      n_clusters == selected_n_clusters
    ) |>
    dplyr::arrange(dplyr::desc(modularity)) |>
    dplyr::slice(1)
  
  best_cluster <- best_run$cluster[[1]]
  best_membership <- best_run$membership[[1]]
  
  membership_df <- data.frame(
    node = names(best_membership),
    cluster = as.integer(best_membership)
  )
  
  output <- list(
    graph = graph,
    runs = leiden_runs,
    summary = leiden_summary,
    stability = leiden_stability,
    eval = leiden_eval,
    candidate_solutions = candidate_solutions,
    selected_resolution = selected_resolution,
    selected_n_clusters = selected_n_clusters,
    best_run = best_run,
    best_cluster = best_cluster,
    best_membership = best_membership,
    membership_df = membership_df
  )
  
  return(output)
}



# ------------------------------------------------------------
# Funktion: Leiden-Diagnoseplots als 3 Panels
# ------------------------------------------------------------


plot_leiden_diagnostics <- function(
    leiden_runs,
    leiden_eval,
    title_prefix = "Leiden",
    smooth = TRUE
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    install.packages("patchwork")
  }
  
  # benötigte Pakete explizit nutzen
  p1 <- ggplot2::ggplot(
    leiden_runs,
    ggplot2::aes(x = resolution, y = n_clusters)
  ) +
    ggplot2::geom_point(alpha = 0.3) +
    {
      if (smooth) ggplot2::geom_smooth(se = FALSE) else NULL
    } +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Resolution",
      y = "Anzahl Cluster",
      title = "A: Clusterzahl entlang Resolution"
    )
  
  p2 <- ggplot2::ggplot(
    leiden_runs,
    ggplot2::aes(x = factor(n_clusters), y = modularity)
  ) +
    ggplot2::geom_boxplot() +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Anzahl Cluster",
      y = "Modularity",
      title = "B: Modularity je Clusterzahl"
    )
  
  p3 <- ggplot2::ggplot(
    leiden_eval,
    ggplot2::aes(
      x = n_clusters,
      y = mean_modularity,
      color = mean_nmi
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Anzahl Cluster",
      y = "Mittlere Modularity",
      color = "Mean NMI",
      title = "C: Güte und Stabilität"
    )
  
  combined_plot <- (p1 / p2 / p3) +
    patchwork::plot_annotation(
      title = title_prefix
    )
  
  return(combined_plot)
}


################################################################################
########### Funktion zum Assessment der wichtigen Variablen für die Cluster-Zuordnung

# ============================================================
# 1. Build cluster interpretation data
# ============================================================

build_cluster_interpretation_data <- function(
    actor_data,
    node_cluster_data,
    module_col_actor = "fk_modul_id",
    actor_col_actor = "Akteur",
    node_col = "name",
    cluster_col = "cluster",
    feature_tables = list(),
    feature_join_col = "fk_modul_id",
    module_name = "modul",
    cluster_prefix = "Cluster_",
    remove_duplicates = TRUE
) {
  
  stopifnot(is.data.frame(actor_data))
  stopifnot(is.data.frame(node_cluster_data))
  
  cluster_ids <- sort(unique(node_cluster_data[[cluster_col]]))
  
  module_cluster_list <- lapply(cluster_ids, function(cl) {
    
    actors_in_cluster <- node_cluster_data[[node_col]][
      node_cluster_data[[cluster_col]] == cl
    ]
    
    modules <- unique(actor_data[[module_col_actor]][
      actor_data[[actor_col_actor]] %in% actors_in_cluster
    ])
    
    data.frame(
      modul = modules,
      clusterinfo = paste0(cluster_prefix, cl),
      stringsAsFactors = FALSE
    )
  })
  
  module_cluster_df <- dplyr::bind_rows(module_cluster_list)
  
  # Achtung: Module können in mehreren Clustern landen,
  # wenn Akteure aus verschiedenen Clustern am selben Modul beteiligt sind.
  module_cluster_df <- dplyr::distinct(module_cluster_df)
  
  feature_data <- dplyr::distinct(
    module_cluster_df[, "modul", drop = FALSE]
  )
  
  for (ft in feature_tables) {
    feature_data <- dplyr::left_join(
      feature_data,
      ft,
      by = stats::setNames(feature_join_col, module_name)
    )
  }
  
  if (remove_duplicates) {
    feature_data <- dplyr::distinct(feature_data)
  }
  
  long_data <- feature_data |>
    tidyr::pivot_longer(
      cols = -dplyr::all_of(module_name),
      names_to = "feature_type",
      values_to = "feature_value"
    ) |>
    dplyr::filter(!is.na(feature_value)) |>
    dplyr::mutate(
      feature_name = paste0(feature_type, "_", feature_value),
      binary = 1
    ) |>
    dplyr::select(
      dplyr::all_of(module_name),
      feature_name,
      binary
    )
  
  wide_data <- tidyr::pivot_wider(
    long_data,
    id_cols = dplyr::all_of(module_name),
    names_from = feature_name,
    values_from = binary,
    values_fn = unique,
    values_fill = 0
  )
  
  wide_data <- as.data.frame(wide_data)
  
  wide_data <- dplyr::left_join(
    wide_data,
    module_cluster_df,
    by = stats::setNames("modul", module_name)
  )
  
  # Random Forest mag saubere Spaltennamen
  old_names <- names(wide_data)
  new_names <- make.names(old_names, unique = TRUE)
  names(wide_data) <- new_names
  
  names_lookup <- data.frame(
    original_name = old_names,
    clean_name = new_names,
    stringsAsFactors = FALSE
  )
  
  list(
    data = wide_data,
    module_cluster_df = module_cluster_df,
    names_lookup = names_lookup
  )
}

####### Cluster interpretieren
# ============================================================
# 2. Assess cluster importance using Random Forest
# ============================================================

assess_cluster_importance_rf <- function(
    data,
    target_cluster,
    cluster_col = "clusterinfo",
    module_col = "modul",
    threshold = 0.75,
    ntree = 500,
    train_fraction = 0.8,
    seed = 123,
    plot = TRUE,
    name_lookup = NULL,
    actor_cluster_data = NULL,
    actor_node_col = "node",
    actor_cluster_col = "cluster",
    verbose = TRUE
) {
  
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package `randomForest` is required.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required.")
  }
  
  data <- as.data.frame(data)
  
  if (!cluster_col %in% names(data)) stop("`cluster_col` not found in `data`.")
  if (!module_col %in% names(data)) stop("`module_col` not found in `data`.")
  if (!target_cluster %in% data[[cluster_col]]) stop("`target_cluster` not found in `cluster_col`.")
  
  data$is_target_cluster <- as.factor(data[[cluster_col]] == target_cluster)
  
  model_data <- data |>
    dplyr::select(-dplyr::all_of(cluster_col))
  
  set.seed(seed)
  
  train_indices <- sample(
    seq_len(nrow(model_data)),
    size = floor(train_fraction * nrow(model_data))
  )
  
  train_data <- model_data[train_indices, ] |>
    dplyr::select(-dplyr::all_of(module_col))
  
  for (i in seq_along(train_data)) {
    if (!is.factor(train_data[[i]])) {
      train_data[[i]] <- as.factor(train_data[[i]])
    }
  }
  
  rf_model <- randomForest::randomForest(
    is_target_cluster ~ .,
    data = train_data,
    ntree = ntree
  )
  
  importance_data <- randomForest::importance(rf_model)
  
  if (is.matrix(importance_data) || is.data.frame(importance_data)) {
    importance_vec <- importance_data[, 1]
  } else {
    importance_vec <- importance_data
  }
  
  importance_vec <- sort(importance_vec, decreasing = TRUE)
  cumulative_sum <- cumsum(importance_vec)
  
  index_threshold <- which(
    cumulative_sum >= threshold * sum(importance_vec)
  )[1]
  
  top_variables <- importance_vec[seq_len(index_threshold)]
  
  result <- data.frame(
    variable = names(top_variables),
    importance = as.numeric(top_variables),
    cumulative_importance = cumsum(top_variables) / sum(importance_vec),
    stringsAsFactors = FALSE
  )
  
  if (!is.null(name_lookup)) {
    
    required_lookup_cols <- c("original_name", "clean_name")
    
    if (!all(required_lookup_cols %in% names(name_lookup))) {
      stop("`name_lookup` must contain columns `original_name` and `clean_name`.")
    }
    
    result$pretty_name <- name_lookup$original_name[
      match(result$variable, name_lookup$clean_name)
    ]
    
    result$pretty_name[is.na(result$pretty_name)] <-
      result$variable[is.na(result$pretty_name)]
    
  } else {
    result$pretty_name <- result$variable
  }
  
  top_features_text <- paste(result$pretty_name, collapse = ", ")
  
  # ------------------------------------------------------------
  # Akteure im Cluster
  # ------------------------------------------------------------
  
  actors_in_cluster <- NULL
  actors_text <- NULL
  
  if (!is.null(actor_cluster_data)) {
    
    actor_cluster_data <- as.data.frame(actor_cluster_data)
    
    if (!actor_node_col %in% names(actor_cluster_data)) {
      stop("`actor_node_col` not found in `actor_cluster_data`.")
    }
    
    if (!actor_cluster_col %in% names(actor_cluster_data)) {
      stop("`actor_cluster_col` not found in `actor_cluster_data`.")
    }
    
    actors_in_cluster <- actor_cluster_data[[actor_node_col]][
      actor_cluster_data[[actor_cluster_col]] == target_cluster
    ]
    
    actors_in_cluster <- sort(unique(actors_in_cluster))
    
    actors_text <- paste(actors_in_cluster, collapse = ", ")
  }
  
  if (verbose) {
    message(
      paste0(
        "Most important variables for ",
        target_cluster,
        ": ",
        top_features_text
      )
    )
    
    if (!is.null(actors_text)) {
      message(
        paste0(
          "Actors in ",
          target_cluster,
          ": ",
          actors_text
        )
      )
    }
  }
  
  if (plot) {
    graphics::barplot(
      top_variables,
      names.arg = result$pretty_name,
      las = 2,
      cex.names = 0.7,
      main = paste("Variable importance -", target_cluster)
    )
  }
  
  output <- list(
    target_cluster = target_cluster,
    model = rf_model,
    importance = result,
    full_importance = importance_vec,
    top_features = result$pretty_name,
    top_features_text = top_features_text,
    actors = actors_in_cluster,
    actors_text = actors_text,
    threshold = threshold,
    ntree = ntree,
    train_fraction = train_fraction
  )
  
  return(output)
}