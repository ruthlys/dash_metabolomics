# load necessary packages
library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dash)
library(dashDaq)
library(dashTable)
library(plotly)
library(heatmaply)
library(dplyr)
library(base64)
library(data.table)
library(vegan)
library(RColorBrewer)
library(colorspace)

####################################################################################################

#################################### APP START $ LOAD FILES ########################################

app <- Dash$new()
# Initiate application

# load files
read_csv <- function(x){
  df <- data.frame(fread(x, header=T))
}

# load files
pca <- read.table("/Users/ruthschmidt/Dropbox/Work/INRS/Data/NI_experiment/Files/VOC/pca.df.txt", sep="\t", header=T)
pca$Date <- as.factor(pca$Date)
df <- read.csv("/Users/ruthschmidt/Dropbox/Work/Plotly/test_files/fc_voc.csv", sep=",", header=T)

# heatmap for ANOVA 0.01 and identified compounds with AMDIS
data.norm <- data.frame(fread("/Users/ruthschmidt/Dropbox/Work/INRS/Data/NI_experiment/Files/VOC/data_normalized.csv", header=T))
# Get anova file 
aov <- data.frame(fread("/Users/ruthschmidt/Dropbox/Work/INRS/Data/NI_experiment/Files/VOC/anova_0.01_cpds_id_final.csv", sep=",", header=T))
# merge and subset data.norm based on aov 
merge.df <- merge(aov, data.norm, by ="ID")
merge.df <- merge.df[-c(1,3:11)]

# make compounds row names
merge.df <- data.frame(merge.df[,-1], row.names=merge.df[,1])
merge.df<-= as.data.frame(t(merge.df))

# add treatment and date from mapping file 
treat <- data.frame(fread("/Users/ruthschmidt/Dropbox/Work/INRS/Data/NI_experiment/Files/VOC/mapping_file.tsv", sep="\t", header=T))
treat <- treat[order(treat$ID), ]
treat <- data.frame(treat[,-1], row.names=treat[,1])
# combine treat and merge.df and set factors 
heat.df <- cbind.data.frame(treat, merge.df)
heat.df$Treatment <- as.factor(heat.df$Treatment)
heat.df$Date <- as.factor(heat.df$Date)
# sort data based on treatment and date
heat.df <- heat.df[order(heat.df$Treatment, heat.df$Date), ]
# gather values to plot in boxplot
box.df <- heat.df %>% tidyr::gather(Compound, value, -Treatment, -Date)

# Define the number of colors 
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Paired"))(nb.cols)

m <- list(
  l = 0,
  r = 0,
  b = 0,
  t = 40
  # pad = 4
)

legend <- list(
  xanchor = "center",
  yanchor = "top",
  # y = "-0.3",
  x = "0.8"
)

create3dScatter <- function(pca_file){
  plot_ly(data = pca_file, x = ~x, y = ~y, z = ~z, 
          color = ~Treatment, colors = "Paired", mode = "markers", type = "scatter3d") %>%
    layout(title = "PCoA of Volatile Organic Compounds (VOCs)", 
           # height=600,width=650,
      scene = list(xaxis = list(title = 'PC1 (37.1%)'),
                        yaxis = list(title = 'PC2 (13.3%)'),
                        zaxis = list(title = 'PC3 (5.7%)'))) %>% layout(margin = m, legend = legend)
  
}

m2 <- list(
  l = 0,
  r = 0,
  b = 20,
  t = 30
  # pad = 4
)

createHeatmap <- function(heat_file){ 
  heatmaply(x = heat_file, 
                         xlab = "Compound", 
                         # ylab = "Sample ID",
                         method = "euclidean",
                         seriate = "mean",
                         row_dend_left = F,
                         plot_method = "plotly",
                         fontsize_col = 10,
                         fontsize_row = 10,
                         key.title='Concentration',
                         showticklabels = c(T,T)) %>% layout(title = "Heatmap of VOCs clustered by Date and Treatment", margin = m2, 
                                                             height=600) %>%
    colorbar(tickfont = list(size = 12), titlefont = list(size = 12), which = 1) %>%
    colorbar(tickfont = list(size = 12), titlefont = list(size = 12), which = 2)
}

t <- list(
  family = "Open Sans",
  size = 12,
  color = "#444444",
  m = list(
    t = "25px"))

dashVolcano <- function (dataframe,  p = "P", effect_size = "EFFECTSIZE", snp = "Compound", 
                                        xlabel = NULL, ylabel = "-log10(p)", point_size = 8,
                                        col = c("#1F78B4"), effect_size_line = c(-1, 1), effect_size_line_color = c("#252525"), effect_size_line_width = 0.5, 
                                        genomewideline_value = -log10(1e-04),
                                        genomewideline_color = "grey", genomewideline_width = 0.5, 
                                        highlight = NULL, highlight_color = c("#1F78B4"), ...) 
{
  manhattanly::volcanoly(x = dataframe, 
                         p = p,
                         snp = snp,
                         effect_size = effect_size,
                         # title = title,
                         xlab = xlabel,
                         ylab = ylabel,
                         point_size = point_size,
                         col = col,
                         effect_size_line = effect_size_line,
                         effect_size_line_color = effect_size_line_color,
                         effect_size_line_width = effect_size_line_width,
                         genomewideline = genomewideline_value,
                         genomewideline_color = genomewideline_color,
                         genomewideline_width = genomewideline_width,
                         highlight = highlight,
                         highlight_color = highlight_color,
                         ...) %>% layout(title = "Volcano Plot of VOCs", font = t)
}

createBox <- function(table_file){
  # create graph parameters 
  plot_ly(
    data = table_file, 
    x = ~Compound, y = ~value,
    color = ~Date, colors = "Paired",
    type = "box", alpha = 0.7) %>% 
    layout(height=500,
      boxmode = "group", 
      title = "Selected VOC across Date", 
      xaxis = list(title = 'Compound', tickangle = 45),
      yaxis = list(title = 'log median norm'))
}

####################################################################################################

#################################### CREATE LAYOUT #################################################

app$layout(
  htmlDiv(
    list(
      # create and style title
      htmlDiv(
        list(
          htmlH2(
            id = "title",
                 'Interactive Metabolomics Viewer', 
                 style = list(
                   textAlign = 'left'
                 )
               ),
          htmlH5(
            'Statistical & functional analysis of metabolomics data', 
            style = list(
              textAlign = 'left'
            )
          )
        ), 
      className = "banner"
      ),
      
      
      
      # create top container and graph child, set style 
      htmlDiv(
        list(
          htmlDiv(
            id = "top-graphs",
            style = list(
              display = "flex", 
              borderRadius = '5px',
              justifyContent = "space-evenly", 
              width = "100%",
              # height = "50px",
              margin = "0 auto"
            ),
            children = list(
              # create tabs
              htmlDiv(
                id = "control-tabs",
                className = "control-tabs",
                children = list(
                  dccTabs(
                    id="tabs", 
                    value = "what-is",
                    children = list(
                      dccTab(
                        label = "About",
                        value = "what-is",
                        children = htmlDiv(
                          className = "control-tab",
                          children = list(
                            htmlH4(
                              className = "what-is", 
                              children = "What is PCoA?"
                            ),
                            dccMarkdown(
                              'Principal Coordinates Analysis (PCoA, = Multidimensional scaling, MDS) 
                      is a method to explore and to visualize similarities or dissimilarities of data. 
                      It starts with a distance matrix annd assigns for each item a 
                      location in a low-dimensional space (here 3D scatter plot).'
                            ),
                            htmlH4(
                              className = "what-is", 
                              children = "What is Volcano Plot?"
                            ),
                            dccMarkdown(
                              'A scatterplot that shows statistical significance (P value) versus magnitude of change (fold change) in this case 
                      of volatile organic compounds (VOCs).'
                            )
                          )
                        )
                      ),
                      
                      dccTab(
                        label = "Parameters",
                        value = "parameters",
                        children = htmlDiv(
                          className = "control-tab",
                          children = list(
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select Date"
                                ),
                                dccDropdown(
                                  id = "d_date-box",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px", 
                                               width = "100%"),
                                  # style = list(width = "45%"),
                                  options = lapply(list("20190723", "20190905"),
                                                   function(x){
                                                     list(label = x, value = x)
                                                   }
                                  ),
                                  value = c('20190723','20190905'),
                                  multi = TRUE
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select Treatment"
                                ),
                                dccDropdown(
                                  id = "d_treatment",
                                  className = "dropdowns",
                                  style = list(marginRight = "10px", width = "100%"),
                                  # style = list(width = "45%"),
                                  options = lapply(list("Control", "NI"),
                                                   function(x){
                                                     list(label = x, value = x)
                                                   }
                                  ),
                                  value = c('Control', 'NI'),
                                  multi = TRUE
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-settings-name",
                                  children = "Settings for Volcano Plot:"
                                ))),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select effect size"
                                ),
                                dccRangeSlider(
                                  id = 'volcanoplot-input',
                                  min = -3,
                                  max = 3,
                                  step = 0.5,
                                  marks = setNames(
                                    lapply(-3:3,
                                           function(i){
                                             list(label = as.character(i))
                                           }),
                                    -3:3
                                  ),
                                  value = c(-0.5, 1)
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                htmlDiv(
                                  className = "app-controls-name",
                                  children = "Select threshold"
                                ),
                                dccSlider(
                                  id = "vp-genomic-line-val",
                                  value = 2,
                                  max = 4,
                                  min = 1,
                                  step = 0.1
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                daqLEDDisplay(
                                  label = "VOCs down",
                                  id = "vp-upper-left",
                                  size = 10,
                                  color = "#19D3F3"
                                )
                              )
                            ),
                            
                            htmlDiv(
                              className = "app-controls-block",
                              children = list(
                                daqLEDDisplay(
                                  label = "VOCs Up",
                                  id = "vp-upper-right",
                                  size = 10,
                                  color = "#19D3F3"
                                )
                              )
                            )
                          )
                        )
                      )
                    )
                  )
                )
              ),
              
              htmlDiv(
                id = "left-top-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "3d-pca",
                        figure = create3dScatter(pca),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container top left
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "2rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "38%", 
                  float = "none", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              ),
              
              htmlDiv(
                id = "right-top-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "volcano-graph",
                        figure = dashVolcano(dataframe = df),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container top right
                style = list(
                  marginTop = "10px",
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "2rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "35%",
                  float = "none",
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              )
            )
          )
        )
      ),
      
      
      # create top container and graph child, set style 
      htmlDiv(
        list(
          htmlDiv(
            id = "bottom-graphs",
            style = list(
              display = "flex",
              borderRadius = '5px',
              justifyContent = "space-evenly", 
              width = "100%",
              margin = "0 auto"
            ),
            children = list(
              htmlDiv(
                id = "left-bottom-graph",
                className = "container",
                list(
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccGraph(
                        id = "heatmap",
                        figure = createHeatmap(heat.df),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container bottom left
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "3rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "48%", 
                  float = "none", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              ),

              htmlDiv(
                id = "right-bottom-graph",
                className = "container",
                list(
                  htmlDiv(
                          className = "app-controls-name",
                          children = "Select Compounds"
                        ),
                  htmlDiv(
                    style = list(width = "100%"),
                    list(
                      dccDropdown(
                        id = "d-cpds-box",
                        className = "dropdowns",
                        style = list(marginTop = "5px", marginBottom = "10px", 
                                     width = "80%"),
                        options = lapply(list("1058_1,8-Octanediol", "1377_alpha-Pinene", "1058_1,8-Octanediol",	"107_1-Hexene, 4-methyl-",	"1072_Unknown",	"1099_Toluene",	"1377_alpha-Pinene",
                                              "279_Glycolaldehyde dimethyl acetal",	"325_Unknown",	"364_Isopropyl acetate",	
                                              "38_Methylacetate", "386_Hexanoic acid methyl ester",	"483_Carbamic acid methyl-, ethylester","	5_1-Hexene, 4-methyl-",	"602_Unknown",	"638_Unknown"),
                                         function(x){
                                           list(label = x, value = x)
                                         }
                        ),
                        value = c('1058_1,8-Octanediol', '1377_alpha-Pinene'),
                        multi = TRUE
                      ),
                      dccGraph(
                        id = "box-graph",
                        figure = createBox(box.df),
                        style = list(width = '100%')
                      )
                    )
                  )
                ),
                # style container top left
                style = list(
                  marginTop = "10px", 
                  marginBottom = "10px",
                  marginLeft = 0,
                  marginRight = 0,
                  paddingTop = "3rem",
                  paddingBottom = "2rem",
                  borderRadius = '5px',
                  width = "48%", 
                  boxSizing = "border-box",
                  boxShadow = '2px 2px 1px #f2f2f2'
                )
              )
            )
          )
        )
      )
    )
  )
)


####################################################################################################

#################################### CALLBACKS #####################################################

app$callback(
  output=list(id="3d-pca", property="figure"),
  list(input(id='d_date-box', property='value'),
       input(id='d_treatment', property='value')),
  function(selected_date, selected_treatment){
    filtered_variable <- pca %>% dplyr::filter(Date %in% c(selected_date), Treatment %in% c(selected_treatment))
    create3dScatter(filtered_variable)
  }
)

app$callback(
  output=list(id="heatmap", property="figure"),
  list(input(id='d_date-box', property='value'),
       input(id='d_treatment', property='value')),
  function(selected_date, selected_treatment){
    filtered_variable <- heat.df %>% dplyr::filter(Date %in% c(selected_date), Treatment %in% c(selected_treatment))
    createHeatmap(filtered_variable)
  }
)

app$callback(
  output = list(id = 'volcano-graph', property = 'figure'),
  params = list(input(id = 'volcanoplot-input', property = 'value'),
                input(id = "vp-genomic-line-val", property = "value")),
  function(effects, genomic_line) {
    dashVolcano(
      dataframe = df,
      genomewideline_value = as.numeric(genomic_line),
      effect_size_line = unlist(effects),
    )
  }
)

app$callback(
  output=list(id="box-graph", property="figure"),
  list(input(id='d_date-box', property='value'),
       input(id='d_treatment', property='value'),
       input(id="d-cpds-box", property="value")),
  function(selected_date, selected_treatment, selected_compound){
    filtered_variable <- box.df %>% dplyr::filter(Date %in% c(selected_date), Treatment %in% c(selected_treatment), 
                                                  Compound %in% c(selected_compound))
    createBox(filtered_variable)
  }
)

app$callback(
  output("vp-upper-right", "value"),
  list(
    input("volcano-graph", "figure"),
    input("vp-genomic-line-val", "value"),
    state("volcanoplot-input", "value")
  ),
  function(fig, thresh, bounds){
    u_lim <- bounds[[2]]
    number = 0
    if (length(fig[["data"]]) > 1){
      x <- unlist(fig[["data"]][[1]]["x"])
      y <- unlist(fig[["data"]][[1]]["y"])
      number <- sum((x > u_lim) & (y > thresh))
    }
    number
  }
)


app$callback(
  output("vp-upper-left", "value"),
  list(
    input("volcano-graph", "figure"),
    input("vp-genomic-line-val", "value"),
    state("volcanoplot-input", "value")
  ),
  function(fig, thresh, bounds){
    u_lim <- bounds[[1]]
    number = 0
    if (length(fig[["data"]]) > 1){
      x <- unlist(fig[["data"]][[1]]["x"])
      y <- unlist(fig[["data"]][[1]]["y"])
      number <- sum((x < u_lim) & (y > thresh))
    }
    number
  }
)


app$run_server(debug=F, threaded=T, showcase = T)

