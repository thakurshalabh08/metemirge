library(ggplot2)
library(RColorBrewer)
library(paletteer)
library(tibble)
library(dplyr)

#### function for paletteer_d reworked from https://github.com/EmilHvitfeldt/paletteer/issues/13 ######

p_d <- function(package, palette, n, direction = 1, type = c("discrete", 
    "continuous")) {
  if (abs(direction) != 1) {
      stop("direction must be 1 or -1")
  }
  #package <- rlang::quo_name(rlang::enquo(package))
  #palette <- rlang::quo_name(rlang::enquo(palette))
  package <- match.arg(package, unique(paletteer::palettes_d_names$package))
  type <- match.arg(type)
  pal <- paletteer::palettes_d[[c(package, palette)]]
  if (is.null(pal)) 
      stop("Palette not found. Make sure the palette name are spelled correct.")
  if (missing(n)) {
      n <- length(pal)
  }
  if (type == "discrete" && n > length(pal)) {
      stop(paste("Number of requested colors greater than this palette can offer which is ", 
          length(pal), ".", sep = ""))
  }
  out <- switch(type, continuous = (grDevices::colorRampPalette(pal))(n), 
      discrete = pal[1:n])
  if (direction == -1) {
      rev(out)
  }
  else {
      out
  }
}

#### function for paletteer_dynamic reworked from https://github.com/EmilHvitfeldt/paletteer/blob/master/R/paletteer_dynamic.R

p_dynamic <- function (package, palette, n, direction = 1) {

  if (abs(direction) != 1) {
    stop("direction must be 1 or -1")
  }


  #package <- rlang::quo_name(rlang::enquo(package))
  #palette <- rlang::quo_name(rlang::enquo(palette))
  if (missing(n)) {
    stop("n not found. Please supply the number of colors you want returned.")
  }

  pal <- paletteer::palettes_dynamic[[c(package, palette)]]
  if (is.null(pal))
    stop("Palette not found. Make sure both package and palette name are spelled correct.")

  if (n > length(pal)) {
    stop(paste("Number of requested colors greater than this palette can offer which is ",
               length(pal), ".", sep = ""))
  }

  if (direction == -1) {
    rev(pal[[n]])
  } else {
    pal[[n]]
  }
}

#### Read taxonomic abundance plot input file from EMIRGE pipeline ####

args=commandArgs(trailingOnly=TRUE)

input_dir<-args[1]
output_file<-args[2]

list_infiles<-list.files(input_dir,recursive=TRUE,full.names=TRUE,pattern="*_abundance_plot_data\\.tsv")

message("Reading abundance results for plotting")
print(list_infiles)

data <- lapply(list_infiles,read.delim) %>% bind_rows

### order data frame with columns sample name  & percent abundance for plotting
plot_data<-data[order(data$SAMPLE_NAME,-data$PERCENT_ABUNDANCE,decreasing=TRUE),]

### read name of the subject species in taxonomy variable #####
message("Reading Taxonomy")
taxonomy<-plot_data$SUBJECT_SPECIES

### parse taxonomy to get genus name and create a data.frame ####
message("Parsing taxonomy to get genus names")
taxa.df=NULL
taxa_list=list()
index_taxa=1

	for(taxa in taxonomy){
		split_taxa<-strsplit(taxa," ")
		taxa.df<-data.frame(taxa,split_taxa[[1]][[1]])
		colnames(taxa.df)<-c("SUBJECT_SPECIES","GENUS")
		taxa_list[[index_taxa]]<-taxa.df
		index_taxa=index_taxa + 1
	}

all_taxa.df<-as.data.frame(do.call("rbind",taxa_list))
colnames(all_taxa.df)<-c("SUBJECT_SPECIES","GENUS")

#### group taxa names by their genus #####
message("Grouping species names by genus")
taxa_group<-all_taxa.df %>% group_by(GENUS) %>% summarize(taxa_cluster=toString(SUBJECT_SPECIES))
species_cluster<-taxa_group$taxa_cluster

#### get color palettee from RColorBrewer ####
#list_palette<-RColorBrewer::brewer.pal.info()
#sequential_palette<-list_palette %>% rownames_to_column("palette_name") %>% filter(list_palette$category=="seq")
#palette_names<-sequential_palette$palette_name
list_palette<-palettes_dynamic_names
sequential_palette<-list_palette %>% filter(list_palette$package=="cartography" && list_palette$type=="sequential")
sequential_palette<-sequential_palette[order(sequential_palette$package,decreasing=FALSE),]
#palette_names<-sequential_palette$palette
#palette_name_list<-as.list(palette_names)
#package_names<-sequential_palette$package
#package_name_list<-as.list(package_names)
#max_color_value<-sequential_palette$length
#max_color_list<-as.list(max_color_value)

#print(palette_name_list)

### get color palettee for genus level from iwanthue function from hues package, will try other palettee options here too ####
message("Assigning color codes to each taxa based on genus")
taxa_color.df=NULL
taxa_color_list=list()
index_taxa_color=1
single_taxa=character(length=0)

for(taxa_row in 1:nrow(taxa_group)){

    species_taxa=taxa_group[taxa_row,2]
    genus_taxa=taxa_group[taxa_row,1]
    species_taxa<-gsub(", ",",",species_taxa)
    split_taxa<-strsplit(species_taxa,",")
    species_taxa2<-unlist(split_taxa)
    ncols<-length(unique(species_taxa2))

 if(ncols==1){
       single_taxa<-append(single_taxa,unique(species_taxa2))
    }else{

	package_name=NULL
	pal_name=NULL
	max_color=NULL

	for(pal_index in 1:nrow(list_palette)){

		 max_color=list_palette[pal_index,3]

                 if(is.null(max_color)){
                    next
                 }

		if(ncols<=max_color){
			package=list_palette[pal_index,1]
			pal_name=list_palette[pal_index,2]
			list_palette<-list_palette[-c(pal_index),]
			break
                }
        }

	#pal_name=palette_name_list[[1]]
	#print(pal_name)
	getPalette = colorRampPalette(p_dynamic(package,pal_name,ncols))
	pal=getPalette(ncols)
	names(pal)=unique(species_taxa2)
	taxa_color.df<-data.frame(names(pal),pal)
	colnames(taxa_color.df)<-c("SUBJECT_SPECIES","COLOR_CODE")
	sub_taxa.df<-all_taxa.df %>% filter(all_taxa.df$SUBJECT_SPECIES %in% species_taxa2)
        taxa_color.df<-merge(taxa_color.df,unique(sub_taxa.df),by=("SUBJECT_SPECIES"),all.x=TRUE,all.y=TRUE)
        colnames(taxa_color.df)<-c("SUBJECT_SPECIES","COLOR_CODE","GENUS")
	taxa_color_list[[index_taxa_color]]<-taxa_color.df
    	index_taxa_color=index_taxa_color+1
    }
}

    ##### add singleton taxa #####
    message("Adding single taxa")
    ncols=length(single_taxa)
    getPalette = colorRampPalette(p_d("ggsci","default_igv",ncols))
    pal=getPalette(ncols)
    names(pal)=unique(single_taxa)
    taxa_color.df<-data.frame(unique(single_taxa),pal)
    colnames(taxa_color.df)<-c("SUBJECT_SPECIES","COLOR_CODE")
    sub_taxa.df<-all_taxa.df %>% filter(all_taxa.df$SUBJECT_SPECIES %in% single_taxa)
    taxa_color.df<-merge(taxa_color.df,unique(sub_taxa.df),by=("SUBJECT_SPECIES"),all.x=TRUE,all.y=TRUE)
    colnames(taxa_color.df)<-c("SUBJECT_SPECIES","COLOR_CODE","GENUS")
    taxa_color_list[[index_taxa_color]]<-taxa_color.df
    

all_taxa_color.df<-as.data.frame(do.call("rbind",taxa_color_list))
colnames(all_taxa_color.df)<-c("SUBJECT_SPECIES","COLOR_CODE","GENUS")

#### merge color code information to main plot data file data.frame by suject species names into new data.frame ####
plot_data2<-merge(plot_data,all_taxa_color.df,by=("SUBJECT_SPECIES"),all.x=TRUE,all.y=TRUE)
#plot_data2<-plot_data2[order(plot_data2$SAMPLE_NAME,-plot_data2$PERCENT_ABUNDANCE),]

### get color code for taxa ####
taxa_color<-as.character(plot_data2$COLOR_CODE)
taxonomy2<-plot_data2$SUBJECT_SPECIES
names(taxa_color)=taxonomy2

### ggplot ###
message("Plotting Relative Abundance for all samples")
message("Writing plot to: ",output_file)
plot<-ggplot(data=plot_data2,aes(x=factor(SAMPLE_NAME),y=plot_data2$PERCENT_ABUNDANCE)) + geom_bar(stat='identity',aes(fill=taxonomy2)) + scale_fill_manual(breaks=plot_data2$SUBJECT_SPECIES, values=taxa_color)
final_plot<-plot + theme_bw() + theme(legend.position = "bottom",legend.direction = "vertical",legend.title.align = 0.5) + ylim(0, 101) + xlab("SAMPLE_NAME") + ylab("RELATIVE ABUNDANCE")

pdf(output_file, width=17,height=12)
print(final_plot)
dev.off()


