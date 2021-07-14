utils::globalVariables(c(
    "peakScore", "geneScore", "x", "y", "label",
    "time_point", "promoter_feature"
))
utils::globalVariables(c("cor", "pvalue"))

#' plot_peakGeneCor
#'
#' @importFrom dplyr filter select as_tibble mutate inner_join
#' @importFrom dplyr group_by summarise left_join
#' @importFrom ggplot2 aes
#' @importFrom stats median
#'
#' @param mmAnnoCor the annotated GRange object from peakGeneCor
#' @param select_gene a gene_id which you want to show
#' @param sigShow one of 'pvalue' 'padj' 'qvalue'
#' @param addLine whether add cor line
#' @param addFullInfo whether add full feature info on plot
#'
#' @return ggplot2 object
#' @export
#'
#' @examples
#'
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     data("RNA_normCount")
#'     data("ATAC_normCount")
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'     peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)[1:100]
#'     mmAnno <- mm_geneScan(peak_GR, Txdb)
#'
#'     ATAC_colData <- data.frame(
#'         row.names = colnames(ATAC_normCount),
#'         type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
#'     )
#'
#'     integrate_replicates(ATAC_normCount, ATAC_colData) -> ATAC_normCount_merge
#'     RNA_colData <- data.frame(
#'         row.names = colnames(RNA_normCount),
#'         type = gsub("_R[0-9]", "", colnames(RNA_normCount))
#'     )
#'     integrate_replicates(RNA_normCount, RNA_colData) -> RNA_normCount_merge
#'     mmAnnoCor <- peakGeneCor(
#'         mmAnno = mmAnno,
#'         peakScoreMt = ATAC_normCount_merge,
#'         geneScoreMt = RNA_normCount_merge,
#'         parallel = FALSE
#'     )
#'
#'     plot_peakGeneCor(mmAnnoCor, select_gene = "AT5G01010")
#'
#' }
plot_peakGeneCor <- function(mmAnnoCor,
                             select_gene,
                             addLine = TRUE,
                             addFullInfo = TRUE,
                             sigShow = c("pvalue", "padj", "qvalue")) {

    sigShow <- match.arg(sigShow, c("pvalue", "padj", "qvalue"))
    peakScoreMt <- metadata(mmAnnoCor)$peakScoreMt
    geneScoreMt <- metadata(mmAnnoCor)$geneScoreMt

    if (!select_gene %in% rownames(geneScoreMt)) {
        stop("the gene you select is not in your geneScoreMt",
            call. = FALSE
        )
    }
    geneExpr <- geneScoreMt[select_gene, ]

    select_df <- mmAnnoCor %>%
        data.frame() %>%
        as_tibble() %>%
        filter(gene_id == select_gene)
    feature_N <- nrow(select_df)

    peak_gene_score <- peakScoreMt[select_df$feature_id, , drop = FALSE] %>%
        as_tibble(rownames = "feature_id") %>%
        tidyr::pivot_longer(
            cols = -c("feature_id"),
            names_to = "time_point",
            values_to = "peakScore"
        ) %>%
        mutate(
            time_point = factor(time_point,
                levels = unique(colnames(peakScoreMt))
            ),
            geneScore = rep(geneExpr, feature_N)
        )

    plot_data <- inner_join(
        peak_gene_score,
        select_df
    )

    p <- ggplot2::ggplot(data = plot_data, aes(
        x = peakScore,
        y = geneScore,
    )) +
        ggplot2::geom_point(alpha = 0.8) +
        ggplot2::facet_wrap(~feature_id, scales = "free") +
        ggplot2::theme_bw()

    if (metadata(mmAnnoCor)$mmCor_mode == "enhancerPromoterCor") {
        select_promoter <- unique(select_df$promoter_feature)
        p <- p +
            ggplot2::labs(
                x = "enhancerScore",
                y = "promoterScore"
            ) +
            ggplot2::ggtitle(paste0(
                "gene: ", select_gene, "\n",
                "promoter: ", select_promoter
            ))
    } else {
        p <- p + ggplot2::ggtitle(select_gene)
    }

    if (addLine) {
        p <- p +
            ggplot2::stat_smooth(
                geom = "line", method = "lm",
                aes(group = feature_id),
                color = "black"
            )
    }

    if (addFullInfo) {
        pos_Info <- plot_data %>%
            group_by(feature_id) %>%
            summarise(
                x = median(peakScore),
                y = median(geneScore)
            ) %>%
            left_join(select_df) %>%
            mutate(label = paste0(
                "dist_to_TSS = ", distanceToTSS,
                "\ncor_value = ", round(cor, digits = 2),
                "\n", sigShow, " = ", format(!!sym(sigShow),
                    scientific = TRUE,
                    digits = 2
                )
            ))

        p <- p +
            ggrepel::geom_text_repel(aes(x = x, y = y, label = label),
                data = pos_Info,
                alpha = 0.8
            )
    }


    return(p)
}


#' shinyParse_peakGeneCor
#'
#' @import shiny
#'
#' @param mmAnnoCor the annotated GRange object from peakGeneCor
#'
#' @return shiny mode
#' @export
#'
#' @examples
#' if (require(TxDb.Athaliana.BioMart.plantsmart28)) {
#'     data("RNA_normCount")
#'     data("ATAC_normCount")
#'     Txdb <- TxDb.Athaliana.BioMart.plantsmart28
#'     seqlevels(Txdb) <- paste0("Chr", c(1:5, "M", "C"))
#'     peak_path <- system.file("extdata", "ATAC.bed.gz", package = "FindIT2")
#'     peak_GR <- loadPeakFile(peak_path)[1:100]
#'     mmAnno <- mm_geneScan(peak_GR, Txdb)
#'
#'     ATAC_colData <- data.frame(
#'         row.names = colnames(ATAC_normCount),
#'         type = gsub("_R[0-9]", "", colnames(ATAC_normCount))
#'     )
#'
#'     ATAC_normCount_merge <- integrate_replicates(ATAC_normCount, ATAC_colData)
#'     RNA_colData <- data.frame(
#'         row.names = colnames(RNA_normCount),
#'         type = gsub("_R[0-9]", "", colnames(RNA_normCount))
#'     )
#'     RNA_normCount_merge <- integrate_replicates(RNA_normCount, RNA_colData)
#'     mmAnnoCor <- peakGeneCor(
#'         mmAnno = mmAnno,
#'         peakScoreMt = ATAC_normCount_merge,
#'         geneScoreMt = RNA_normCount_merge,
#'         parallel = FALSE
#'     )
#'
#'     # shinyParse_peakGeneCor(mmAnnoCor)
#'
#' }
shinyParse_peakGeneCor <- function(mmAnnoCor) {
    all_genes <- unique(mmAnnoCor$gene_id)

    mmAnnoCor$distanceToTSS <- as.integer(mmAnnoCor$distanceToTSS)

    ui <- shiny::fluidPage(
        fluidRow(
            column(
                2,
                suppressWarnings(selectizeInput("gene",
                    "Gene",
                    choices = all_genes
                )),
            ),
            column(
                2,
                numericInput("cor", "Cor_filter",
                    value = 0, min = 0, max = 1,
                    step = 0.01
                )
            ),
            column(
                2,
                numericInput("pvalue", "pvalue_filter",
                    value = 1, min = 0, max = 1,
                    step = 1e-50
                )
            )
        ),
        fluidRow(
            column(
                4,
                tableOutput("GRange"),
                tableOutput("geneScore"),
                tableOutput("peakScore")
            )
        ),
        fluidRow(
            column(
                12,
                actionButton("click", "plot"),
                actionButton("reset", "Clear"),
                plotOutput("plot_peakGene")
            )
        )
    )


    peakScoreMt <- metadata(mmAnnoCor)$peakScoreMt %>%
        tibble::as_tibble(rownames = "feature_id")
    geneScoreMt <- metadata(mmAnnoCor)$geneScoreMt %>%
        tibble::as_tibble(rownames = "gene_id")

    server <- function(input, output, session) {
        data <- reactive(subset(
            mmAnnoCor,
            gene_id %in% input$gene &
                abs(cor) > input$cor &
                pvalue < input$pvalue
        ))

        output$GRange <- renderTable(data.frame(data()))
        output$peakScore <- renderTable(
            dplyr::filter(peakScoreMt, feature_id %in% data()$feature_id)
        )
        output$geneScore <- renderTable(geneScoreMt %>%
            dplyr::filter(gene_id == input$gene))

        # https://shiny.rstudio.com/articles/action-buttons.html
        v <- reactiveValues(data = NULL)
        observeEvent(input$click, {
            v$data <- plot_peakGeneCor(
                mmAnnoCor = data(),
                select_gene = input$gene
            )
        })

        observeEvent(input$reset, {
            v$data <- NULL
        })


        output$plot_peakGene <- renderPlot(
            {
                if (is.null(v$data)) {
                    return()
                }
                v$data
            },
            res = 96
        )
    }

    shinyApp(ui, server)
}
