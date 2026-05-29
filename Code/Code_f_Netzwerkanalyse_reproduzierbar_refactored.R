################################################################################
# Reproduzierbare Netzwerkanalyse der Akteursstruktur im bundesweiten Monitoring
#
# Zweck:
#   Dieses Skript erzeugt aus einer Modul-Akteur-Tabelle ein Kooperationsnetzwerk,
#   identifiziert robuste Leiden-Cluster, interpretiert diese fachlich über
#   Moduleigenschaften und erstellt eine publikationsfähige statische Abbildung.
#
# Voraussetzungen:
#   Die folgenden Objekte müssen vor Ausführung geladen sein:
#   - Akteure:      Tabelle mit mindestens fk_modul_id, Akteur, Akteur_lang
#   - alleModule:   Vektor der zu berücksichtigenden Modul-IDs
#   - Kopfdaten:    Metadaten der Module, u.a. modul_id, Progrshort
#   - BiolFokus:    Modultabelle mit fk_modul_id, Sphaere, Fokus
#   - ArtengrDat:   Modultabelle mit fk_modul_id, Artengruppe, BV_Art
#   - HabitateDat:  Modultabelle mit fk_modul_id, habitat, BV_habitat
#
# Wichtige ausgelagerte Funktionen:
#   - make_contingency()
#   - run_leiden_cluster_workflow()
#   - plot_leiden_diagnostics()
#   - build_cluster_interpretation_data()
#   - assess_cluster_importance_rf()
#
# Hinweis:
#   Das Skript trennt Analyse- und Darstellungsentscheidungen. Manuelle
#   Koordinatenanpassungen oder das Entfernen von Ausreißern aus Clusterhüllen
#   wirken nur auf die Abbildung, nicht auf Cluster- oder Zentralitätsberechnungen.
################################################################################

# ==============================================================================
# 0. Setup
# ============================================================================== 

required_packages <- c(
  "tidyverse", "igraph", "ggraph", "tidygraph", "ggforce", "ggnewscale",
  "patchwork", "cowplot", "networkD3", "openxlsx", "randomForest", "scales",
  "grid", "plotly", "svglite"
)

invisible(lapply(required_packages, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package not installed: ", pkg, call. = FALSE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

# Pfade bei Bedarf anpassen -----------------------------------------------------
helper_file <- "Helper_functions_netzwerkanalyse.R"
output_dir  <- "outputs/netzwerkanalyse"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(helper_file)

# Reproduzierbarkeit der Layouts und Stichproben
set.seed(123)

# Prüfen, ob die zentralen Eingabeobjekte vorhanden sind ------------------------
required_objects <- c(
  "Akteure", "alleModule", "Kopfdaten", "BiolFokus", "ArtengrDat", "HabitateDat"
)

missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects) > 0) {
  stop(
    "Folgende Eingabeobjekte fehlen im Workspace: ",
    paste(missing_objects, collapse = ", "),
    call. = FALSE
  )
}

# ==============================================================================
# 1. Daten vorbereiten
# ============================================================================== 

# Es werden nur Module berücksichtigt, die auch in der Flowchart-Auswahl enthalten
# sind. Jede Zeile beschreibt hier: Organisation X ist an Modul Y beteiligt.
network <- Akteure |>
  dplyr::filter(fk_modul_id %in% alleModule) |>
  dplyr::select(Project = fk_modul_id, Org = Akteur) |>
  dplyr::distinct(Project, Org, .keep_all = TRUE)

# Kurzer Überblick über berücksichtigte Programme/Module.
selected_programmes <- unique(Kopfdaten$Progrshort[Kopfdaten$modul_id %in% alleModule])
print(selected_programmes)

# ---------------------------------------------------------------------------- 
# 1.1 Harmonisierung von Akteurskürzeln
# ---------------------------------------------------------------------------- 
# Ziel: fachlich identische oder für die Abbildung zusammenzufassende Einheiten
# werden unter einem gemeinsamen Kürzel geführt. Das reduziert künstliche
# Fragmentierung im Netzwerk.

network <- network |>
  dplyr::mutate(
    Org = dplyr::case_when(
      Org %in% c("TI-SF", "TI-WO", "TI-OF") ~ "TI",
      Org %in% c("Jagdverbände") ~ "DJV",
      Org %in% c("Senckenberg-Wild", "Senckenberg-Naturkunde", "Senckenberg-Meer") ~ "Senckenberg",
      Org %in% c("WI-SSG", "WI-CRG") ~ "WI",
      Org %in% c("ICES-WGBEAM", "ICES-WGBIFS") ~ "ICES",
      TRUE ~ Org
    )
  )

# Entfernen unscharfer Sammelkategorien oder fachlich für diese Auswertung wenig
# hilfreicher Einträge. Die Bereinigung ist eine analytische Entscheidung und
# sollte bei Veröffentlichung dokumentiert werden.
excluded_orgs <- c(
  "Bingo", "IngoBrandt", "LUPUS", "Nationalparkämter", "Länder", "BIMA",
  "OEKOKART", "PAN", "Vogelschutzwarten", "naturschutzfonds-BB",
  "Ehrenamt-MV", "Feuchtgebiete", "NABU-AGKranich", "WI", "IWSG"
)

network <- network |>
  dplyr::filter(!Org %in% excluded_orgs)

# ==============================================================================
# 2. Akteursgruppen für Visualisierung definieren
# ============================================================================== 

Akteure_Uebersicht <- data.frame(
  Akteur_kurz = unique(network$Org),
  Akteur_lang = Akteure$Akteur_lang[match(unique(network$Org), Akteure$Akteur)],
  stringsAsFactors = FALSE
)

# Achtung: Diese Zuordnung ist positionsbasiert und daher empfindlich gegenüber
# Änderungen in der Reihenfolge von unique(network$Org). Für langfristige
# Reproduzierbarkeit wäre eine externe Mapping-Tabelle vorzuziehen.
akteur_kategorien <- c(
  "Bund","Bund","Wissenschaft","Länder","Bund",NA,NA,"Fachverbände",NA,
  "Wissenschaft","Wissenschaft","Länder","Wissenschaft","Länder",
  "Fachverbände","Wissenschaft","Fachverbände","Wissenschaft",
  "Fachverbände","Fachverbände","Länder","Fachverbände","Wissenschaft",
  "Länder","Länder","Länder","Länder","Länder","Länder","Länder",
  "Bund","Bund","Wissenschaft",NA,"Länder","Bund","Wissenschaft",
  "Länder","Wissenschaft","Länder","Länder","Fachverbände","Länder",
  "Fachverbände","Fachverbände","Länder","Länder","Länder",
  "Wissenschaft","Wissenschaft","Wissenschaft","Länder","Bund",
  "Fachverbände","Länder","Fachverbände","Wissenschaft","Länder",
  "Wissenschaft","Fachverbände","Länder","Länder","Länder","Länder",
  "Länder","Länder","Länder","Fachverbände","Länder","Wissenschaft",
  "Fachverbände","Fachverbände","Wissenschaft","Fachverbände",
  "Fachverbände","Fachverbände","Fachverbände","Fachverbände",
  "Länder","Länder","Länder","Fachverbände","Länder","Länder",
  "Länder","Länder","Länder","Länder","Länder","Länder","Länder",
  "Länder","Länder","Länder","Länder","Länder","Länder","Länder",
  "Länder","Länder","Länder","Länder","Länder","Bund","Länder",
  "Länder","Länder","Länder","Länder","Länder","Länder","Länder",
  "Länder","Länder","Länder","Länder","Länder","Länder"
)

if (length(akteur_kategorien) != nrow(Akteure_Uebersicht)) {
  warning(
    "Die Länge des manuell definierten Kategorie-Vektors passt nicht zur Anzahl ",
    "der Akteure. Bitte Akteure_Uebersicht prüfen. Fehlende Kategorien werden NA."
  )
  Akteure_Uebersicht$Kategorie <- NA_character_
  Akteure_Uebersicht$Kategorie[seq_along(akteur_kategorien)] <- akteur_kategorien
} else {
  Akteure_Uebersicht$Kategorie <- akteur_kategorien
}

# OSPAR/HELCOM u.a. werden hier der Gruppe Bund zugeschlagen. Das ist fachlich
# eine vereinfachende Visualisierungsentscheidung.
Akteure_Uebersicht$Kategorie[is.na(Akteure_Uebersicht$Kategorie)] <- "Bund"

farben <- c(
  "Bund" = "#fcd500",
  "Länder" = "#0fa446",
  "Wissenschaft" = "#a2d109",
  "Fachverbände" = "#95d2ec"
)

Akteure_Uebersicht$Farbe <- farben[Akteure_Uebersicht$Kategorie]

# ==============================================================================
# 3. Netzwerk aus gemeinsamer Modulbeteiligung aufbauen
# ============================================================================== 

# make_contingency() kommt aus Helper_functions_netzwerkanalyse.R.
# Sie erzeugt alle Paarungen von Organisationen, die gemeinsam an einem Modul
# beteiligt sind.
kontingenztab <- make_contingency(
  dataset = network,
  Project = "Project",
  Organization = "Org"
)

# Kantengewicht: Anzahl gemeinsamer Module je Organisationspaar.
edge_data <- kontingenztab |>
  dplyr::group_by(V1, V2) |>
  dplyr::summarise(n_Projs = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(V1, V2) |>
  dplyr::filter(V1 != V2) |>
  dplyr::mutate(weight = scales::rescale(n_Projs, to = c(1, 10)))

network_full <- igraph::graph_from_data_frame(
  d = edge_data,
  directed = FALSE
)

# Für die Abbildung und Clusteranalyse wird die größte zusammenhängende Komponente
# verwendet, damit isolierte Kleinstkomponenten die Struktur nicht dominieren.
components_full <- igraph::components(network_full)
giant_component <- which.max(components_full$csize)
nodes_to_keep <- igraph::V(network_full)[components_full$membership == giant_component]
network_main <- igraph::induced_subgraph(network_full, vids = nodes_to_keep)

# ==============================================================================
# 4. Leiden-Clusteranalyse und Diagnose
# ============================================================================== 

# Die optimale globale Clusterlösung wird über verschiedene Resolution-Werte und
# wiederholte Läufe anhand von Modularity, Stabilität (NMI) und Clusterzahl geprüft.
global_leiden <- run_leiden_cluster_workflow(
  graph = network_main,
  resolution_values = seq(0.2, 3, by = 0.05),
  n_runs = 100,
  min_clusters = 5,
  seed = 123,
  use_weights = TRUE
)

p_leiden_global <- plot_leiden_diagnostics(
  leiden_runs = global_leiden$runs,
  leiden_eval = global_leiden$eval,
  title_prefix = "Leiden-Diagnostik: Gesamtnetzwerk"
)

print(p_leiden_global)

membership_df <- global_leiden$membership_df

# Clusterattribute an igraph hängen.
igraph::V(network_main)$cluster <- membership_df$cluster[
  match(igraph::V(network_main)$name, membership_df$node)
]

# ==============================================================================
# 5. Fachliche Interpretation der globalen Cluster
# ============================================================================== 

# Datensatz für Random-Forest-Interpretation der Cluster vorbereiten.
cluster_interpretation <- build_cluster_interpretation_data(
  actor_data = Akteure,
  node_cluster_data = membership_df,
  module_col_actor = "fk_modul_id",
  actor_col_actor = "Akteur",
  node_col = "node",
  cluster_col = "cluster",
  feature_tables = list(
    BiolFokus[c("fk_modul_id", "Sphaere", "Fokus")],
    ArtengrDat[c("fk_modul_id", "Artengruppe", "BV_Art")],
    HabitateDat[c("fk_modul_id", "habitat", "BV_habitat")]
  ),
  feature_join_col = "fk_modul_id",
  module_name = "modul",
  cluster_prefix = "Cluster_"
)

wide_data <- cluster_interpretation$data

cluster_results <- purrr::map(
  unique(wide_data$clusterinfo),
  \(cl) assess_cluster_importance_rf(
    data = wide_data,
    target_cluster = cl,
    cluster_col = "clusterinfo",
    module_col = "modul",
    threshold = 0.75,
    ntree = 500,
    train_fraction = 0.8,
    seed = 123,
    plot = TRUE,
    name_lookup = cluster_interpretation$names_lookup,
    actor_cluster_data = membership_df |>
      dplyr::mutate(cluster = paste0("Cluster_", cluster)),
    actor_node_col = "node",
    actor_cluster_col = "cluster",
    verbose = TRUE
  )
)

names(cluster_results) <- unique(wide_data$clusterinfo)

cluster_importance_table <- purrr::map_dfr(
  cluster_results,
  \(x) x$importance |> dplyr::mutate(cluster = x$target_cluster)
) |>
  dplyr::arrange(cluster, dplyr::desc(importance))

cluster_actor_table <- purrr::map_dfr(
  cluster_results,
  \(x) data.frame(cluster = x$target_cluster, actor = x$actors)
)

# Optional exportieren.
openxlsx::write.xlsx(
  list(
    variable_importance = cluster_importance_table,
    actors_by_cluster = cluster_actor_table
  ),
  file = file.path(output_dir, "cluster_interpretation_global.xlsx"),
  overwrite = TRUE
)

# ==============================================================================
# 6. Subclusteranalyse innerhalb des marinen/Nationalpark-Clusters
# ============================================================================== 

# Fachliche Beobachtung: Der globale Cluster 1 enthält marine Akteure und
# Nationalpark-Akteure. Deshalb wird dieser Cluster separat als Subgraph geprüft.
marine_parent_cluster <- 1

marine_nodes <- membership_df |>
  dplyr::filter(cluster == marine_parent_cluster) |>
  dplyr::pull(node)

marine_subgraph <- igraph::induced_subgraph(network_main, vids = marine_nodes)

marine_leiden <- run_leiden_cluster_workflow(
  graph = marine_subgraph,
  resolution_values = seq(0.2, 3, by = 0.05),
  n_runs = 100,
  min_clusters = 2,
  seed = 123,
  use_weights = TRUE
)

p_leiden_marine <- plot_leiden_diagnostics(
  leiden_runs = marine_leiden$runs,
  leiden_eval = marine_leiden$eval,
  title_prefix = "Leiden-Diagnostik: mariner Subcluster"
)

print(p_leiden_marine)

marine_membership_df <- marine_leiden$membership_df |>
  dplyr::rename(subcluster = cluster)

marine_subcluster_interpretation <- build_cluster_interpretation_data(
  actor_data = Akteure,
  node_cluster_data = marine_membership_df,
  module_col_actor = "fk_modul_id",
  actor_col_actor = "Akteur",
  node_col = "node",
  cluster_col = "subcluster",
  feature_tables = list(
    BiolFokus[c("fk_modul_id", "Sphaere", "Fokus")],
    ArtengrDat[c("fk_modul_id", "Artengruppe", "BV_Art")],
    HabitateDat[c("fk_modul_id", "habitat", "BV_habitat")]
  ),
  feature_join_col = "fk_modul_id",
  cluster_prefix = "Marine_Subcluster_"
)

marine_wide_data <- marine_subcluster_interpretation$data

marine_subcluster_results <- purrr::map(
  unique(marine_wide_data$clusterinfo),
  \(cl) assess_cluster_importance_rf(
    data = marine_wide_data,
    target_cluster = cl,
    threshold = 0.75,
    name_lookup = marine_subcluster_interpretation$names_lookup,
    actor_cluster_data = marine_membership_df |>
      dplyr::mutate(subcluster = paste0("Marine_Subcluster_", subcluster)),
    actor_node_col = "node",
    actor_cluster_col = "subcluster",
    plot = TRUE,
    verbose = TRUE
  )
)

names(marine_subcluster_results) <- unique(marine_wide_data$clusterinfo)

marine_subcluster_importance_table <- purrr::map_dfr(
  marine_subcluster_results,
  \(x) x$importance |> dplyr::mutate(cluster = x$target_cluster)
)

marine_subcluster_actor_table <- purrr::map_dfr(
  marine_subcluster_results,
  \(x) data.frame(cluster = x$target_cluster, actor = x$actors)
)

openxlsx::write.xlsx(
  list(
    variable_importance = marine_subcluster_importance_table,
    actors_by_subcluster = marine_subcluster_actor_table
  ),
  file = file.path(output_dir, "cluster_interpretation_marine_subclusters.xlsx"),
  overwrite = TRUE
)

# ==============================================================================
# 7. Finale Clusterzuordnung für die Abbildung
# ============================================================================== 

# Für die Visualisierung wird der marine/Nationalpark-Cluster in zwei Cluster
# aufgeteilt. Damit wird aus der hierarchischen Struktur eine flache Darstellung:
#   1 = Marines Monitoring
#   6 = Monitoring der Nationalparke
# Die anderen globalen Cluster bleiben unverändert.

membership_plot <- membership_df |>
  dplyr::mutate(cluster_plot = cluster)

marine_subcluster_1_actors <- marine_subcluster_results$Marine_Subcluster_1$actors
marine_subcluster_2_actors <- marine_subcluster_results$Marine_Subcluster_2$actors

# Falls die Interpretation der beiden Subcluster fachlich andersherum ausfällt,
# können die folgenden zwei Zeilen einfach vertauscht werden.
membership_plot$cluster_plot[
  membership_plot$cluster == marine_parent_cluster &
    membership_plot$node %in% marine_subcluster_1_actors
] <- 1

membership_plot$cluster_plot[
  membership_plot$cluster == marine_parent_cluster &
    membership_plot$node %in% marine_subcluster_2_actors
] <- 6

cluster_names <- c(
  "1" = "Marines Monitoring",
  "2" = "Forstliches Monitoring",
  "3" = "Naturschutzmonitoring,\nbehördliche Akteure (z.B. FFH, WRRL)",
  "4" = "Naturschutzmonitoring,\nehrenamtliche Akteure (z.B. Vogelmonitoring)",
  "5" = "edaphisches Monitoring,\nwissenschaftliche Akteure",
  "6" = "Monitoring der Nationalparke"
)

cluster_hull_colors <- c(
  "1" = "#66C2A5",
  "2" = "#FFD92F",
  "3" = "#FC8D62",
  "4" = "#E78AC3",
  "5" = "#A6D854",
  "6" = "#8DA0CB"
)

# ==============================================================================
# 8. Layout, Zentralität und Plotdaten vorbereiten
# ============================================================================== 

# Zentralität wird auf dem vollständigen finalen Netzwerk berechnet, nicht auf
# manuell verschobenen Plot-Koordinaten.
igraph::V(network_main)$cluster_plot <- membership_plot$cluster_plot[
  match(igraph::V(network_main)$name, membership_plot$node)
]

igraph::V(network_main)$Kategorie <- Akteure_Uebersicht$Kategorie[
  match(igraph::V(network_main)$name, Akteure_Uebersicht$Akteur_kurz)
]

# Gewichtete Betweenness: starke Kooperationsbeziehungen werden als kurze Wege
# interpretiert, daher 1 / weight.
igraph::V(network_main)$betweenness <- igraph::betweenness(
  network_main,
  normalized = TRUE,
  weights = 1 / igraph::E(network_main)$weight
)

igraph::V(network_main)$betweenness_scaled <- scales::rescale(
  igraph::V(network_main)$betweenness,
  to = c(2.5, 10)
)

graph_tbl <- tidygraph::as_tbl_graph(network_main)

layout_df <- ggraph::create_layout(
  graph_tbl,
  layout = "nicely"
)

# Die layout_tbl_graph-Klasse muss erhalten bleiben, damit ggraph() später noch
# Kanten und Knoten korrekt kennt. Deshalb werden Plotattribute per Basiszugriff
# ergänzt statt das Layout in ein normales data.frame umzuwandeln.
layout_df$cluster <- factor(layout_df$cluster_plot, levels = names(cluster_names))

# Organisationen, die wegen der konvexen Hüllen optional nicht in die Hüllen
# eingehen sollen. Die Punkte bleiben sichtbar, nur die Hulls ignorieren sie.
outlier_orgs_for_hulls <- c(
  "MEROS", "GCD", "ICES", "Senckenberg-Meer", "ITAW", "Uni-Trier",
  "SBMS-HB", "LfL-BY", "Uni-Hamburg"
)

# Manuelle Koordinatenanpassungen für die finale Abbildung. Diese Änderungen sind
# rein grafisch und verändern keine Netzwerkmetriken.
manual_offsets <- dplyr::tribble(
  ~name,               ~dx,  ~dy,
  "SBMS-HB",          -4.5,  0,
  "MEROS",            -2.2,  1.94,
  "GCD",              -1.3,  1.5,
  "ICES",              0,    1.94,
  "Senckenberg-Meer",  0,    1.72,
  "Uni-Trier",         0,    1.3,
  "LfL-BY",            0.5, -0.3
)

layout_df_plot <- layout_df

for (i in seq_len(nrow(manual_offsets))) {
  idx <- layout_df_plot$name == manual_offsets$name[i]
  layout_df_plot$x[idx] <- layout_df_plot$x[idx] + manual_offsets$dx[i]
  layout_df_plot$y[idx] <- layout_df_plot$y[idx] + manual_offsets$dy[i]
}

# Hüllen auf Basis der finalen Plot-Koordinaten berechnen. Für die Hüllen reicht
# ein normales data.frame; der eigentliche Netzwerkplot nutzt weiterhin das
# layout_tbl_graph-Objekt layout_df_plot.
hull_df_plot <- as.data.frame(layout_df_plot) |>
  dplyr::filter(!is.na(cluster), !name %in% outlier_orgs_for_hulls)

# Labels: pro Cluster die obersten 20 % nach Betweenness plus manuell ergänzte
# Schlüsselorganisationen.
extra_label_nodes <- c("UBA")

label_nodes <- as.data.frame(layout_df_plot) |>
  dplyr::filter(!is.na(cluster)) |>
  dplyr::group_by(cluster) |>
  dplyr::arrange(dplyr::desc(betweenness), .by_group = TRUE) |>
  dplyr::mutate(
    rank = dplyr::row_number(),
    n_cluster = dplyr::n(),
    keep = rank <= ceiling(0.20 * n_cluster)
  ) |>
  dplyr::filter(keep) |>
  dplyr::ungroup() |>
  dplyr::pull(name) |>
  unique()

label_nodes <- unique(c(label_nodes, extra_label_nodes))

layout_df_plot$label <- ifelse(
  layout_df_plot$name %in% label_nodes,
  layout_df_plot$name,
  NA_character_
)

# ==============================================================================
# 9. Ausreißerdiagnose auf Basis der Plot-Koordinaten
# ============================================================================== 

outlier_candidates <- as.data.frame(layout_df_plot) |>
  dplyr::filter(!is.na(cluster)) |>
  dplyr::group_by(cluster) |>
  dplyr::mutate(
    cluster_x = mean(x, na.rm = TRUE),
    cluster_y = mean(y, na.rm = TRUE),
    dist_to_cluster_center = sqrt((x - cluster_x)^2 + (y - cluster_y)^2)
  ) |>
  dplyr::arrange(cluster, dplyr::desc(dist_to_cluster_center)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::ungroup() |>
  dplyr::select(
    cluster, name, Kategorie, x, y, cluster_x, cluster_y,
    betweenness, dist_to_cluster_center
  )

print(outlier_candidates, n = 60)

# ==============================================================================
# 10. Statischer Plot mit separater Clusterlegende
# ============================================================================== 

network_plot_static <- ggraph::ggraph(layout_df_plot) +
  ggforce::geom_mark_hull(
    data = hull_df_plot,
    ggplot2::aes(x = x, y = y, group = cluster, fill = cluster),
    color = "grey45",
    alpha = 0.10,
    concavity = 5,
    expand = grid::unit(5, "mm"),
    radius = grid::unit(3, "mm"),
    linetype = "dashed",
    linewidth = 0.7,
    show.legend = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_hull_colors,
    labels = cluster_names,
    name = "Netzwerkcluster",
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +
  ggraph::geom_edge_link(
    color = "grey75",
    alpha = 0.45,
    width = 0.4
  ) +
  ggraph::geom_node_point(
    ggplot2::aes(
      size = betweenness_scaled,
      fill = Kategorie,
      text = paste0(
        "Organisation: ", name,
        "<br>Cluster: ", cluster,
        "<br>Betweenness: ", round(betweenness, 3),
        "<br>x: ", round(x, 3),
        "<br>y: ", round(y, 3)
      )
    ),
    shape = 21,
    color = "grey25",
    stroke = 0.45,
    alpha = 0.9
  ) +
  ggplot2::scale_fill_manual(
    values = farben,
    name = "Akteursgruppe",
    na.value = "grey80",
    guide = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        shape = 21,
        color = "grey25",
        linetype = 0,
        linewidth = 0,
        alpha = 1,
        size = 4
      ),
      order = 1
    )
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.5, 10),
    name = "Brückenfunktion",
    breaks = c(2.5, 6.25, 10),
    labels = c("gering", "mittel", "stark"),
    guide = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0,
      override.aes = list(
        shape = 21,
        fill = "grey70",
        color = "grey25",
        linetype = 0,
        linewidth = 0,
        alpha = 1
      ),
      order = 2
    )
  ) +
  ggraph::geom_node_text(
    ggplot2::aes(label = label),
    repel = TRUE,
    size = 3,
    color = "grey10",
    max.overlaps = Inf,
    box.padding = grid::unit(0.4, "lines"),
    point.padding = grid::unit(0.3, "lines")
  ) +
  ggplot2::theme_void() +
  ggplot2::labs(
    title = "Akteursnetzwerk im bundesweiten Biodiversitätsmonitoring"
  ) +
  ggplot2::theme(
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "top",
    plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.text = ggplot2::element_text(size = 9)
  )

print(network_plot_static)

# Separater Dummy-Plot nur für die Clusterlegende unten.
cluster_legend_df <- data.frame(
  cluster = factor(names(cluster_names), levels = names(cluster_names)),
  x = 1,
  y = 1
)

cluster_legend_plot <- ggplot2::ggplot(
  cluster_legend_df,
  ggplot2::aes(x = x, y = y, fill = cluster)
) +
  ggplot2::geom_point(
    ggplot2::aes(x = x, y = y, fill = cluster),
    shape = 22,
    size = 5,
    color = "grey45",
    alpha = 0,
    show.legend = TRUE
  ) +
  ggplot2::scale_fill_manual(
    values = cluster_hull_colors,
    labels = cluster_names,
    name = "Netzwerkcluster",
    guide = ggplot2::guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        shape = 22,
        size = 5,
        color = "grey45",
        linetype = "dashed",
        linewidth = 0.8,
        alpha = 0.35
      )
    )
  ) +
  ggplot2::theme_void() +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_text(size = 10, face = "bold"),
    legend.text = ggplot2::element_text(size = 9)
  )

legend_components <- cowplot::get_plot_component(
  cluster_legend_plot,
  "guide-box",
  return_all = TRUE
)

cluster_legend <- legend_components[[which.max(
  vapply(legend_components, function(x) {
    as.numeric(grid::convertWidth(grid::grobWidth(x), "cm"))
  }, numeric(1))
)]]

network_plot_final <- network_plot_static / cluster_legend +
  patchwork::plot_layout(heights = c(1, 0.16))

print(network_plot_final)

# ==============================================================================
# 11. Optional: interaktive Plotly-Version zur Prüfung von Punktnamen
# ============================================================================== 

plotly_plot <- plotly::ggplotly(network_plot_static, tooltip = "text")
# plotly_plot

# ==============================================================================
# 12. Export
# ============================================================================== 

ggplot2::ggsave(
  filename = file.path(output_dir, "Netzwerkplot_Cluster.png"),
  plot = network_plot_final,
  width = 3200,
  height = 2400,
  units = "px",
  dpi = 300
)

ggplot2::ggsave(
  filename = file.path(output_dir, "Netzwerkplot_Cluster.svg"),
  plot = network_plot_final,
  width = 3200,
  height = 2400,
  units = "px",
  device = svglite::svglite
)

# Ende -------------------------------------------------------------------------
