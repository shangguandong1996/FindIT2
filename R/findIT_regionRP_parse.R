utils::globalVariables(c(".", "percent", "merge_result"))

#' shinyParse_findIT_regionRP
#'
#' if you have regionRP result from calcRP_region, and TF_result from
#' findIT_regionRP. you can just c(regionRP, TF_result) to merge two result. Then
#' use shinyParse_findIT_regionRP to parse the result
#'
#' @param merge_result the merge result from c(regionRP, TF_result)
#' @param mode one of 'gene' or 'TF'
#'
#' @return shiny website
#' @export
#'
shinyParse_findIT_regionRP <- function(merge_result,
                                       mode = "gene"){

    mode <- match.arg(mode, c("TF", "gene"))

    sumRP_data <- SummarizedExperiment::assays(merge_result)$sumRP
    fullRP_data <- SummarizedExperiment::assays(merge_result)$fullRP

    fullRP_GR <- SummarizedExperiment::rowRanges(
        MultiAssayExperiment::experiments(merge_result)[[1]]
    )
    fullRP_GR$distanceToTSS <- as.integer(fullRP_GR$distanceToTSS)

    percent_data <- metadata(merge_result)$percent_df


    percent_data_wider <- metadata(merge_result)$percent_df %>%
        tidyr::pivot_wider(names_from = sample,
                           values_from = percent)


    hits_df <- metadata(merge_result)$hits_df

    TF_percent <- SummarizedExperiment::assays(merge_result)$TF_percent
    TF_pvalue <- SummarizedExperiment::assays(merge_result)$TF_pvalue
    TF_pvalue_rank <- apply(TF_pvalue, 2, rank)

    all_genes <- unique(percent_data$gene_id)
    all_TF <- rownames(TF_percent)
    all_samples <- colnames(TF_percent)

    percent_data$sample <- factor(percent_data$sample, levels = colnames(TF_percent))

    if (mode == "gene"){
        payAttention_gene(all_genes = all_genes,
                          all_TF = all_TF,
                          fullRP_GR = fullRP_GR,
                          fullRP_data = fullRP_data,
                          sumRP_data = sumRP_data,
                          percent_data_wider = percent_data_wider,
                          hits_df = hits_df)
    } else {
        payAttention_TF(all_genes = all_genes,
                        all_TF = all_TF,
                        all_samples = all_samples,
                        TF_pvalue = TF_pvalue,
                        TF_pvalue_rank = TF_pvalue_rank,
                        TF_percent = TF_percent,
                        percent_data = percent_data)
    }

}


payAttention_gene <- function(all_genes,
                              all_TF,
                              fullRP_data,
                              fullRP_GR,
                              sumRP_data,
                              percent_data_wider,
                              hits_df){
    ui <- fluidPage(
        fluidRow(
            column(2,selectizeInput("gene","Gene",choices = all_genes)),
            column(2,selectInput("tf", "TF", choices = all_TF))
        ),
        fluidRow(
            column(8,
                   tableOutput("sumRP"),
                   tableOutput("gene_percent"),
                   tableOutput("fullRP"),
                   tableOutput("GRange")
            )
        )
    )

    server <- function(input, output, session) {
        data <- reactive(subset(fullRP_GR, gene_id == input$gene))

        output$sumRP <- renderTable(sumRP_data[input$gene, , drop = FALSE] %>%
                                        tibble::as_tibble(rownames = "gene_id"))

        output$gene_percent <- renderTable(percent_data_wider %>%
                                               dplyr::filter(
                                                   gene_id == input$gene,
                                                   TF_id == input$tf)
        )



        output$fullRP <- renderTable(cbind(mcols(data())[, c("feature_id",
                                                             "gene_id")],
                                           fullRP_data[names(data()), ]) %>%
                                         tibble::as_tibble() %>%
                                         dplyr::inner_join(
                                             dplyr::filter(hits_df,
                                                           TF_id == input$tf)
                                         ) %>%
                                         dplyr::select(feature_id,
                                                       TF_id,
                                                       dplyr::everything()) %>%
                                         dplyr::select(-gene_id)
        )

        output$GRange <- renderTable(data.frame(data()))

    }

    shinyApp(ui, server)
}

#' @importFrom stats setNames
payAttention_TF <- function(all_genes,
                            all_TF,
                            all_samples,
                            TF_pvalue,
                            TF_pvalue_rank,
                            TF_percent,
                            percent_data
){

    ui <- fluidPage(
        fluidRow(
            column(2,selectInput("tf", "TF", choices = all_TF))
        ),
        fluidRow(
            column(8,
                   tableOutput("TF_pvalue"),
                   tableOutput("TF_rank"),
                   tableOutput("TF_mean_percent"),
                   tableOutput("TF_percent_summary")
            )
        ),
        fluidRow(
            column(12,
                   selectInput("sample", "sample_selct", choices = all_samples,
                               multiple = TRUE),
                   actionButton("click", "plot"),
                   actionButton("reset", "Clear"),
                   plotOutput("percent_distribution"))
        )

    )

    server <- function(input, output, session) {
        output$TF_pvalue <- renderTable(TF_pvalue[input$tf, , drop = FALSE] %>%
                                            tibble::as_tibble(rownames = "TF_id"),
                                        digits = -2)

        output$TF_rank <- renderTable(TF_pvalue_rank[input$tf, , drop = FALSE] %>%
                                          tibble::as_tibble(rownames = "TF_id"))

        output$TF_mean_percent <- renderTable(TF_percent[input$tf, , drop = FALSE] %>%
                                                  tibble::as_tibble(rownames = "TF_id"),
                                              digits = 3)

        output$TF_percent_summary <- renderTable(percent_data %>%
                                                     dplyr::filter(TF_id == input$tf) %>%
                                                     dplyr::group_split(sample) %>%
                                                     purrr::map(function(x) {
                                                         summary_data <- summary(x$percent)
                                                         sample_name <- unique(x$sample)
                                                         summary_name <- names(summary_data)
                                                         matrix(as.numeric(summary_data),nrow = 1) %>%
                                                             as.data.frame() %>%
                                                             setNames(summary_name) %>%
                                                             dplyr::mutate(sample_name = sample_name) %>%
                                                             dplyr::select(sample_name, dplyr::everything())
                                                     }) %>%
                                                     do.call(rbind, .),
                                                 digits = 3)



        # https://shiny.rstudio.com/articles/action-buttons.html
        v <- reactiveValues(data = NULL)
        observeEvent(input$click,{
            v$data <- percent_data %>%
                filter(TF_id == input$tf,
                       sample == input$sample) %>%
                ggplot2::ggplot(aes(x = percent,
                                    color = sample)) +
                ggplot2::geom_histogram(position="identity",
                                        fill = "white",
                                        binwidth = 0.02,
                                        alpha = 0.5,
                                        boundary= 0) +
                ggplot2::theme_bw()
        })

        observeEvent(input$reset, {
            v$data <- NULL
        })

        output$percent_distribution <- renderPlot({
            if (is.null(v$data)) return()
            v$data
        }, res = 96)


    }

    shinyApp(ui, server)

}
