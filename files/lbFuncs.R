## Print infinite rows of a tibble
print.inf <- function(tibble){
  print(tibble, n = Inf)
}

## Get gene symbol from gene ID
idToSym <- function(geneID, genesGR){
  genesGR[match(geneID, genes$gene_id)]$gene_name
}

## Get gene ID from gene symbol
symToId <- function(geneSym, genesGR = genes){
  genesGR[match(geneSym, genes$gene_name)]$gene_id
}

## Common axis labels
log10plab <- expression(paste(-log[10], "(p)"))
log2cpmlab <- expression(paste(log[2], "(CPM)"))

## Remove AsIs class
unAsIs <- function(x){
  if ("AsIs" %in% class(x)) {
    class(x) <- class(x)[-match("AsIs", class(x))]
  }
  x
}

## Simple DT::datatable
## - for small tables needing simple presentation
lbSimpleDT <- function(tbl, cap, rownames = FALSE){
  datatable(
    tbl,
    caption = htmltools::tags$caption(
      cap,
      style = "text-align: left;"
    ),
    options = list(
      dom = "t",  # Available options: tflipr
      columnDefs = list(list(className = 'dt-left', targets = "_all")),  # Align all columns left
      pageLength = -1,  # Show all results
      scrollCollapse = TRUE,  # Removes whitespace when table is small
      scrollY = 350,  # Sets max height of table before scrollbar is added
      scrollX = TRUE,  # Enable horizontal scrolling when too wide
      sScrollX = "100%"  # Stops last column title being cut off
    ),
    rownames = rownames,
    height = "100%",  # Seems to work in tandem with scrollCollapse
    width = "100%",  # Fills width when viewed in RStudio viewer
    class = "stripe hover row-border order-column compact nowrap"
  )
}

## Complex DT::datatable
## - for large tables needing pagination and search
lbComplexDT <- function(tbl, cap, rownames = FALSE){
  datatable(
    tbl,
    caption = htmltools::tags$caption(
      cap,
      style = "text-align: left;"
    ),
    options = list(
      dom = "tflipr",  # Available options: tflipr
      columnDefs = list(list(className = 'dt-left', targets = "_all")),  # Align all columns left
      pageLength = 10,  # Show all results
      scrollCollapse = TRUE,  # Removes whitespace when table is small
      scrollY = "100%",  # Sets max height of table before scrollbar is added
      scrollX = TRUE,  # Enable horizontal scrolling when too wide
      sScrollX = "100%"  # Stops last column title being cut off
    ),
    rownames = rownames,
    height = "100%",  # Seems to work in tandem with scrollCollapse
    width = "100%",  # Fills width when viewed in RStudio viewer
    class = "stripe hover row-border order-column compact nowrap"
  )
}

# Manhattan plot
lbManPlot <- function(tt, highlightChr){
  chrs <- tt$chromosome %>%
    unique() %>%
    mixedsort()
  man <- tt %>%
    mutate(mid = (start + end) / 2) %>%
    dplyr::select(gene_id, gene_name, chromosome, mid, FDR) %>%
    dplyr::filter(!(chromosome %in% c("X", "Y", "MT"))) %>%
    group_by(chromosome) %>%
    summarise(chrLen = max(mid)) %>%
    mutate(chrSt = cumsum(chrLen)-chrLen) %>%
    dplyr::select(-chrLen) %>%
    left_join(tt) %>%
    mutate(mid = (start + end) / 2) %>%
    mutate(midCum = chrSt + mid) %>%
    mutate(
      colour = ifelse(chromosome %in% chrs[c(TRUE, FALSE)], "grey50", "grey80"),
      colour = ifelse(chromosome == highlightChr, "red", colour)
    )
  axis <- man %>%
    group_by(chromosome) %>%
    summarise(center = (max(midCum) + min(midCum)) / 2)
  ggplot(man, aes(x = midCum, y = -log10(PValue))) +
    geom_point(
      data = man %>%
        dplyr::filter(!DE),
      aes(colour = colour),
      shape = 20,
      alpha = 1,
      size = 1
    ) +
    scale_colour_identity() +
    geom_point(
      data = man %>%
        dplyr::filter(DE),
      aes(fill = colour),
      colour = "black",
      shape = 23,
      alpha = 1,
      size = 1,
      stroke = 0.3
    ) +
    scale_fill_identity() +
    scale_x_continuous(label = axis$chromosome, breaks = axis$center) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(
      x = "Chromosome",
      y = expression(paste(-log[10], "(p)"))
    ) +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text = element_text(size = rel(0.8)),
      axis.title = element_text(size = rel(0.8))
    )
}
