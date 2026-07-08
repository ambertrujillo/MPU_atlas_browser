# ============================================================
# CRITICAL: this must be the absolute first code that runs.
# Setting RETICULATE_PYTHON before library(reticulate) is loaded
# prevents reticulate from silently auto-binding to the default
# r-reticulate virtualenv before your explicit path is honored.
# Once ANY reticulate call fires, the Python binding is locked
# for the life of the R session — no code after this point can
# override it. If you ever see the "already been initialized"
# error again, this ordering is the first thing to check.
# ============================================================
library(shiny)
library(ggplot2)
library(ggpubr)
library(anndata)
library(reticulate)


# ============================================================
# Config — pseudobulk CSV lookup
# ============================================================
h5ad_lookup <- list(
  thalamus = list(
    full_atlas      = "data/combined_reclustered_trimmed.h5ad",
    neuron_subtypes = "data/Neu_adata_trimmed.h5ad"
  )
)

adata_atlas_choices <- c("Thalamus" = "thalamus")

adata_subset_choices <- c(
  "Full Atlas"           = "full_atlas",
  "Neuron Subtypes Only" = "neuron_subtypes"
)

# UNCONFIRMED: placeholder filenames — replace with your actual SVG names,
# and confirm these files live in www/ (Shiny only serves static assets
# from www/; a file sitting in data/ will not resolve via tags$img(),
# regardless of a correct-looking relative path).
svg_lookup <- list(
  thalamus = list(
    full_atlas      = "umap_combined_Neu.svg",
    neuron_subtypes = "umap_Neu.svg"
  )
)

file_lookup <- list(
  thalamus = list(
    sampleid          = "data/Thalamus_SampleID.csv",
    sampleid_celltype = "data/Thalamus_generalDisease_celltype.csv"
  )
)

pseudobulk_choices <- c(
  "SampleID"               = "sampleid",
  "SampleID and Cell type" = "sampleid_celltype"
)

xaxis_choices <- c(
  "Disease"             = "general_disease",
  "Disease (split FTD)" = "general_disease_split_FTD"
)

celltype_colors <- c(
  "CP" = "#d5c18e", "BVC" = "#74858c", "OPC" = "#ed7020", "OL" = "#fdb515",
  "Astrocyte" = "#abd158", "PVM" = "#231f20", "Lymphocyte" = "#68151f",
  "Microglia" = "#d24a75", "ExNeu RNF220-" = "#90C6E2", "ExNeu RNF220+" = "#547cb2",
  "InNeu SOX14+" = "#21409a", "InNeu SOX14-" = "#262262"
)
# NOTE: hardcoded, manually-maintained palette. If a new celltype is added
# to the underlying annotation and this vector isn't updated,
# scale_fill_manual() silently drops that celltype's bar rather than erroring.

BASIS <- "X_umap_mnn"
LAYER <- "cpm"

.adata_cache <- new.env()
.adata_cache$adata         <- NULL
.adata_cache$loaded_atlas  <- NULL
.adata_cache$loaded_subset <- NULL
.adata_cache$status        <- "not_loaded"
.adata_cache$error_msg     <- NULL

# ============================================================
# UI
# ============================================================
ui <- fluidPage(
  
  tags$img(src = "card_mpu_logo.svg", height = "100px"),
  
  h4("Pseudobulk expression browser"),
  selectInput("atlas_choice", "Select atlas:",
              choices = c("Thalamus" = "thalamus")),
  selectInput("pseudobulk_level", "Pseudobulked by:",
              choices = pseudobulk_choices),
  actionButton("load_btn", "Load Data"),
  uiOutput("load_status"),
  hr(),
  
  conditionalPanel(
    condition = "input.pseudobulk_level == 'sampleid'",
    selectInput("xaxis_choice", "Plot by:", choices = xaxis_choices),
    selectizeInput("comparison_groups", "Groups to compare (Wilcox test):",
                   choices = NULL, multiple = TRUE,
                   options = list(placeholder = "No comparison"))
  ),
  
  textInput("gene_id", "Enter Gene ID:", value = ""),
  actionButton("plot_btn", "Generate Plot"),
  plotOutput("expr_violin"),
  plotOutput("expr_celltype_bar"),
  
  hr(),
  h4("Feature plot (single-cell UMAP) browser"),
  fluidRow(
    column(6,
           selectInput("adata_atlas_choice", "Select atlas:",
                       choices = adata_atlas_choices),
           selectInput("adata_subset_choice", "AnnData subset:",
                       choices = adata_subset_choices),
           actionButton("load_adata_btn", "Load AnnData (one-time, cached)"),
           uiOutput("adata_status")
    ),
    column(6,
           uiOutput("adata_reference_img")
    )
  ),
  textInput("feature_gene_id", "Enter Gene ID (feature plot):", value = ""),
  actionButton("feature_plot_btn", "Generate Feature Plot"),
  plotOutput("feature_umap")
)

# ============================================================
# Server
# ============================================================
server <- function(input, output, session) {
  
  # Quit the entire R process when the browser tab/session closes.
  # Intended for local, single-user execution (per app.sh). If this app
  # is ever deployed on Shiny Server for multiple concurrent users, this
  # needs to be removed — otherwise one person closing their tab kills
  # the app for everyone else connected.
  session$onSessionEnded(function() {
    stopApp()
  })
  
  # ---- Pseudobulk CSV loading ----
  loaded_df  <- reactiveVal(NULL)
  load_state <- reactiveVal("idle")
  load_label <- reactiveVal(NULL)
  load_error <- reactiveVal(NULL)
  adata_tick <- reactiveVal(0)
  

  displayed_atlas  <- reactiveVal(NULL)
  displayed_subset <- reactiveVal(NULL)
  
  observeEvent(input$load_btn, {
    req(input$atlas_choice, input$pseudobulk_level)
    
    load_state("loading")
    load_error(NULL)
    
    withProgress(message = "Loading data...", value = 0, {
      
      incProgress(0.2, detail = "Resolving file path")
      file_path <- file_lookup[[input$atlas_choice]][[input$pseudobulk_level]]
      
      result <- tryCatch({
        if (is.null(file_path)) {
          stop("No file configured for atlas: ", input$atlas_choice,
               " / level: ", input$pseudobulk_level)
        }
        if (!file.exists(file_path)) {
          stop("File missing from disk: ", file_path)
        }
        
        incProgress(0.5, detail = "Reading CSV")
        df <- read.csv(file_path, stringsAsFactors = FALSE)
        
        incProgress(0.3, detail = "Done")
        df
      }, error = function(e) e)
      
      if (inherits(result, "error")) {
        load_state("error")
        load_error(conditionMessage(result))
        loaded_df(NULL)
      } else {
        loaded_df(result)
        load_label(paste0(input$atlas_choice, " / ",
                          names(pseudobulk_choices)[pseudobulk_choices == input$pseudobulk_level]))
        load_state("success")
      }
    })
  })
  
  output$load_status <- renderUI({
    state <- load_state()
    if (state == "idle") {
      tags$span(style = "color: grey;", "No data loaded yet.")
    } else if (state == "loading") {
      tags$span(style = "color: #cc9900;", "⏳ Loading...")
    } else if (state == "error") {
      tags$span(style = "color: #cc0000;", paste("✗ Load failed:", load_error()))
    } else if (state == "success") {
      tags$span(style = "color: #228822; font-weight: bold;",
                paste("✓ Loaded:", load_label(), "(", nrow(loaded_df()), "rows )"))
    }
  })
  
  observeEvent({load_state(); input$xaxis_choice}, {
    req(load_state() == "success", input$pseudobulk_level == "sampleid", input$xaxis_choice)
    df <- loaded_df()
    
    validate(need(input$xaxis_choice %in% names(df),
                  paste("Column", input$xaxis_choice, "not in loaded data.")))
    
    grp_levels <- sort(unique(as.character(df[[input$xaxis_choice]])))
    
    updateSelectizeInput(session, "comparison_groups",
                         choices = grp_levels, selected = character(0),
                         server = TRUE)
  }, ignoreInit = TRUE)
  
  # ---- Sample-level violin plot ----
  gene_data <- eventReactive(input$plot_btn, {
    req(input$gene_id, input$xaxis_choice)
    
    validate(
      need(load_state() == "success",
           "Load data successfully before generating a plot."),
      need(input$pseudobulk_level == "sampleid",
           "This plot requires 'SampleID' pseudobulk level."),
      need(length(input$comparison_groups) != 1,
           "Select either 0 groups (no comparison) or 2+ groups to compare.")
    )
    
    withProgress(message = "Generating plot...", value = 0, {
      
      incProgress(0.2, detail = "Using loaded data")
      df <- loaded_df()
      
      validate(
        need(input$gene_id %in% df$gene,
             paste("Gene", input$gene_id, "not found in this dataset."))
      )
      
      incProgress(0.3, detail = "Filtering gene")
      out <- subset(df, gene == input$gene_id)
      
      present <- unique(as.character(out[[input$xaxis_choice]]))
      missing_groups <- setdiff(input$comparison_groups, present)
      validate(
        need(length(missing_groups) == 0,
             paste("Selected group(s) not present for this gene:",
                   paste(missing_groups, collapse = ", ")))
      )
      
      attr(out, "gene_label")  <- input$gene_id
      attr(out, "xaxis_col")   <- input$xaxis_choice
      attr(out, "xaxis_label") <- names(xaxis_choices)[xaxis_choices == input$xaxis_choice]
      attr(out, "comparisons") <- if (length(input$comparison_groups) >= 2) {
        combn(input$comparison_groups, 2, simplify = FALSE)
      } else {
        list()
      }
      
      incProgress(0.5, detail = "Rendering")
      out
    })
  })
  
  output$expr_violin <- renderPlot({
    req(input$pseudobulk_level == "sampleid")
    
    df        <- gene_data()
    xcol      <- attr(df, "xaxis_col")
    xlabel    <- attr(df, "xaxis_label")
    genelabel <- attr(df, "gene_label")
    comps     <- attr(df, "comparisons")
    
    group_colors <- c(
      "AD" = "#ff6600", "Control" = "#6666ff",
      "FTD" = "#669900", "FTD-GRN" = "#669900", "sFTD-TDP" = "#99ff66"
    )
    
    p <- ggplot(df, aes(x = .data[[xcol]], y = expression, fill = .data[[xcol]])) +
      geom_violin() +
      geom_point(aes(color = study),
                 position = position_jitterdodge(dodge.width = 0),
                 size = 1.5) +
      stat_summary(fun = "mean", geom = "crossbar", width = 0.25, colour = "red") +
      scale_fill_manual(values = group_colors) +
      theme_bw() +
      labs(y = "Normalized Expression", x = xlabel, title = genelabel) +
      theme(axis.title.x = element_blank())
    
    if (length(comps) > 0) {
      p <- p + stat_compare_means(comparisons = comps, method = "wilcox.test",
                                  label = "p.format")
    }
    
    p
  })
  
  # ---- Celltype-level barplot ----
  # NOTE: this file is pre-aggregated (one row per gene x celltype x disease,
  # no SampleID column) — no donor-level variance is recoverable from this
  # CSV, so no error bars are possible without re-exporting from per-donor
  # data. Flagged, not fixed.
  celltype_data <- eventReactive(input$plot_btn, {
    req(input$gene_id)
    
    validate(
      need(load_state() == "success",
           "Load data successfully before generating a plot."),
      need(input$pseudobulk_level == "sampleid_celltype",
           "This plot requires 'SampleID and Cell type' pseudobulk level.")
    )
    
    withProgress(message = "Generating plot...", value = 0, {
      incProgress(0.2, detail = "Using loaded data")
      df <- loaded_df()
      
      validate(
        need(input$gene_id %in% df$gene,
             paste("Gene", input$gene_id, "not found in this dataset.")),
        need("general_disease" %in% names(df),
             "Expected column 'general_disease' not found in this file."),
        need("celltype" %in% names(df),
             "Expected column 'celltype' not found in this file.")
      )
      
      incProgress(0.4, detail = "Filtering gene")
      out <- subset(df, gene == input$gene_id)
      
      incProgress(0.4, detail = "Rendering")
      attr(out, "gene_label") <- input$gene_id
      out
    })
  })
  
  output$expr_celltype_bar <- renderPlot({
    req(input$pseudobulk_level == "sampleid_celltype")
    
    df        <- celltype_data()
    genelabel <- attr(df, "gene_label")
    
    ggplot(df, aes(x = reorder(celltype, -expression), y = expression, fill = celltype)) +
      geom_col() +
      scale_fill_manual(values = celltype_colors) +
      labs(y = paste0("Normalized Expression (", genelabel, ")")) +
      theme_bw() +
      facet_wrap(~general_disease) +
      theme(legend.title = element_blank(),
            legend.text = element_text(size = 12),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 12),
            axis.title.y = element_text(size = 12))
  })
  
  # ---- AnnData load-once-per-process, cached (single-slot version) ----
  get_adata <- function(requested_atlas, requested_subset) {
    
    same_selection <- identical(.adata_cache$loaded_atlas, requested_atlas) &&
      identical(.adata_cache$loaded_subset, requested_subset)
    
    if (.adata_cache$status == "ready" && same_selection) {
      return(.adata_cache$adata)
    }
    
    if (.adata_cache$status == "loading") {
      while (.adata_cache$status == "loading") Sys.sleep(0.5)
      if (same_selection) return(.adata_cache$adata)
    }
    
    if (.adata_cache$status == "ready" && !same_selection) {
      .adata_cache$adata <- NULL
      gc()
    }
    
    .adata_cache$status <- "loading"
    adata_tick(isolate(adata_tick()) + 1)
    
    file_path <- h5ad_lookup[[requested_atlas]][[requested_subset]]
    
    withProgress(
      message = paste("Loading", names(adata_atlas_choices)[adata_atlas_choices == requested_atlas],
                      "-", names(adata_subset_choices)[adata_subset_choices == requested_subset], "..."),
      value = 0, {
        
        result <- tryCatch({
          incProgress(0.1, detail = "Checking file")
          if (is.null(file_path)) {
            stop("No file configured for atlas '", requested_atlas,
                 "' / subset '", requested_subset, "'")
          }
          if (!file.exists(file_path)) {
            stop("File not found: ", file_path)
          }
          
          incProgress(0.2, detail = "Reading h5ad")
          adata <- read_h5ad(file_path)
          
          incProgress(0.6, detail = "Caching")
          adata
        }, error = function(e) e)
        
        if (inherits(result, "error")) {
          .adata_cache$status        <- "error"
          .adata_cache$error_msg     <- conditionMessage(result)
          .adata_cache$loaded_atlas  <- NULL
          .adata_cache$loaded_subset <- NULL
          message("AnnData load FAILED: ", conditionMessage(result))
          adata_tick(isolate(adata_tick()) + 1)
          return(NULL)
        }
        
        .adata_cache$adata         <- result
        .adata_cache$loaded_atlas  <- requested_atlas
        .adata_cache$loaded_subset <- requested_subset
        .adata_cache$status        <- "ready"
        incProgress(0.1, detail = "Done")
      })
    
    adata_tick(isolate(adata_tick()) + 1)
    .adata_cache$adata
  }
  
  output$adata_status <- renderUI({
    adata_tick()
    atlas_label  <- names(adata_atlas_choices)[adata_atlas_choices == .adata_cache$loaded_atlas]
    subset_label <- names(adata_subset_choices)[adata_subset_choices == .adata_cache$loaded_subset]
    
    switch(.adata_cache$status,
           "not_loaded" = tags$span(style="color:grey;", "AnnData not yet loaded."),
           "loading"    = tags$span(style="color:#cc9900;", "⏳ Loading AnnData..."),
           "ready"      = tags$span(style="color:#228822;font-weight:bold;",
                                    paste("✓ Loaded:", atlas_label, "-", subset_label,
                                          "—", .adata_cache$adata$n_obs, "nuclei")),
           "error"      = tags$span(style="color:#cc0000;", paste("✗ Load failed:", .adata_cache$error_msg))
    )
  })
  
  observeEvent(input$load_adata_btn, {
    # Set immediately on click — this is what makes the reference image
    # swap instantly, decoupled from whether the load actually succeeds.
    displayed_atlas(input$adata_atlas_choice)
    displayed_subset(input$adata_subset_choice)
    
    message(">>> LOAD CLICKED: atlas=", input$adata_atlas_choice,
            " subset=", input$adata_subset_choice, " at ", Sys.time())
    get_adata(input$adata_atlas_choice, input$adata_subset_choice)
  })
  
  output$adata_reference_img <- renderUI({
    req(displayed_atlas(), displayed_subset())
    
    svg_path <- svg_lookup[[displayed_atlas()]][[displayed_subset()]]
    
    validate(
      need(!is.null(svg_path),
           paste("No reference image configured for", displayed_atlas(), "/", displayed_subset()))
    )
    
    tags$img(src = svg_path, style = "max-width: 600px; height: auto;")
  })
  
  feature_data <- eventReactive(input$feature_plot_btn, {
    req(input$feature_gene_id, input$adata_atlas_choice, input$adata_subset_choice)
    
    validate(
      need(.adata_cache$status == "ready",
           "Load an AnnData object first (button above)."),
      need(identical(.adata_cache$loaded_atlas, input$adata_atlas_choice) &&
             identical(.adata_cache$loaded_subset, input$adata_subset_choice),
           "Currently loaded object doesn't match your current selection. Click Load again.")
    )
    
    adata <- get_adata(input$adata_atlas_choice, input$adata_subset_choice)
    
    withProgress(message = "Generating feature plot...", value = 0, {
      incProgress(0.1, detail = "Validating object contents")
      
      validate(
        need(BASIS %in% names(adata$obsm),
             paste("Embedding", BASIS, "not found. Available:", paste(names(adata$obsm), collapse = ", "))),
        need(LAYER %in% adata$layers$keys(),
             paste("Layer", LAYER, "not found. Available:", paste(adata$layers$keys(), collapse = ", "))),
        need(input$feature_gene_id %in% adata$var_names,
             paste("Gene", input$feature_gene_id, "not found in this object."))
      )
      
      incProgress(0.3, detail = "Locating gene and embedding coordinates")
      gene_idx <- which(adata$var_names == input$feature_gene_id)
      
      validate(
        need(length(gene_idx) == 1,
             paste("Gene", input$feature_gene_id, "matched", length(gene_idx),
                   "entries in var_names — expected exactly 1"))
      )
      
      coords <- adata$obsm[[BASIS]]
      
      incProgress(0.4, detail = "Extracting expression values")
      expr <- as.numeric(adata$layers$get(LAYER)[, gene_idx])
      
      incProgress(0.2, detail = "Preparing plot data")
      df <- data.frame(UMAP1 = coords[, 1], UMAP2 = coords[, 2], expression = expr)
      df <- df[order(df$expression), ]
      
      attr(df, "gene_label") <- input$feature_gene_id
      df
    })
  })
  
  output$feature_umap <- renderPlot({
    df <- feature_data()
    genelabel <- attr(df, "gene_label")
    
    ggplot(df, aes(x = UMAP1, y = UMAP2, color = expression)) +
      geom_point(size = 0.5) +
      scale_color_gradientn(colors = c("yellow", "red", "purple", "darkblue")) +
      coord_fixed() +
      theme_void() +
      labs(title = genelabel, color = "Expression")
  })
}

shinyApp(ui, server)