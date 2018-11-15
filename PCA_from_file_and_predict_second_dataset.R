#' PCA and predict second dataset
#'
#' Builds PLS model from training dataset and predicts second dataset
#' Writes out predicted.scores of the second dataset
#' Plots original PLS, projected samples only, and projected samples ontop of original PLS
#'
#' @param file file for X matirx
#' @param sample.names Vector of sample names in X matrix
#' @param sample.type vector of sample groups
#' @param y.response numeric vector of response values in same order as samples
#' @param comp number of components to compute
#' @param scale default=T
#' @param labels label the plot, default = T
#' @param comp.x,comp.y comps to display
#' @param title title of the plot
#' 
#' @param file2 file for test data matrix
#' @param rotate rotates graph along y axis
#' @param sample.names2 Vector of sample names in 2nd dataset, if needed
#' @param sample.type2 Vector of sample types in 2nd dataset, if needed
#' @param train_string string of training data to insert in file name of predicted scores
#' @param test_string string of data being tested to insert in file name of predicted scores
#' @param saveplot whether to save the plot, default is F
#' @param savetype the type of plot to save,options are ".pdf" or ".png"
#' @param w is width of plot to be saved
#' @param h is height of plot to be saved
#' @param legendname is the legend name
#' @param plot_both if true plots both training and test set in color
#' @param colpalette allows you to put in a color palette of form c("#F87660", "#39B600",....etc) to manually assign colors
#' @param shape.palette allows you to put in a shape palette of form c(1, 3,....etc) to manually assign shapes
#' @param varimax If T performs Varimax rotation, 
#' @param varimax.comp # of varimax components, kind of hacky, keep this # the same as # of comps. Will fix later.
#' @param output_folder the folder to output to, default is ./ i.e. current folder
#' @param TCGA predicted files are from TCGA, barcodes separated by periods, so remove normal samples, default is FALSE
#'
#' @importFrom mixOmics pls
#' @export
#'

PCA_from_file_and_predict_second_dataset=function (file, file2, sample.names, sample.type, y.response, 
                                                    sample.names2 = NULL, sample.type2 = NULL, train_string, 
                                                    test_string, title = "PCA", comp.x = "PC1", comp.y = "PC2", 
                                                    comps = 2, labels = F, saveplot = T, savetype = ".png", w = 8, 
                                                    h = 6, legendname = "default", scale = F, center=T, plot_both = T, 
                                                    colpalette = NULL, shape.palette = NULL, ellipses = T, conf = 0.9, 
                                                    varimax = F, varimax.comp = 2, output_folder = "./",TCGA=F,threshold=3,rotate=F,do.legend=T) {
  

  data = read.table(file, sep = "\t", header = T, stringsAsFactors = FALSE,   quote = "")
  data = data[rowSums((data[, -1] == 0)) < ncol(data[-1]), ]
  data2 = read.table(file2, sep = "\t", header = T, stringsAsFactors = FALSE, quote = "")
  
  
  if(TCGA == T){
    temp_name=colnames(data2)[1]
    cancer_samples=which(as.numeric(sapply(colnames(data2)[-1],function(x) strsplit(x,"\\.")[[1]][4])) <= 9) 
    data2=cbind(data2[,1],data2[,-1][,cancer_samples])
    colnames(data2)[1]=temp_name
  }
  data = data[!duplicated(data[, 1]), ]
  data2 = data2[!duplicated(data2[, 1]), ]
  common.genes = intersect_all(data[, 1], data2[, 1])
  data = data[data[, 1] %in% common.genes, ]
  data2 = data2[data2[, 1] %in% common.genes, ]
  data = data[order(data[, 1]), ]
  data2 = data2[order(data2[, 1]), ]
  rownames(data) = make.names(data[, 1], unique = TRUE)
  t.data = data.frame(t(data[, -1]))
  pca<-prcomp(t.data,scale=scale,center=center,rank.=comps);
  x.variates=as.data.frame(pca$x)
 # x.variates=cbind("Sample"=rownames(x.variates),x.variates)
  x.loadings=pca$rotation
  #pca_loadings=cbind("Loading"=rownames(pca_loadings),pca_loadings)#if genenames in rownames
  #x.loadings=cbind("Loading"=data[,1],x.loadings)#if genenames not in rownames
  pca_evalues=pca$sdev
  
  
  
  rownames(data2) = make.names(data2[, 1], unique = TRUE)
  t.data2 = data.frame(t(data2[, -1]))
  
  
  
  rotated.data2 = scale(t.data2, pca$center, pca$scale) %*% as.matrix(x.loadings) 
  
  
  
  
  
  
  
  
  if (varimax == T) {
    rotation = varimax(as.matrix(x.loadings[, c(1:(varimax.comp))]), 
                       normalize = F)
    scores <- as.matrix(x.variates[, 1:(varimax.comp)]) %*% 
      rotation$rotmat
    scores = as.data.frame(scores)
    colnames(scores)[1:ncol(scores)]=colnames(x.variates)[1:ncol(scores)]
    x.variates = scores
    x.loadings=loadings(rotation)[,1:(varimax.comp)]
  }
  
  if(rotate==T){
    temp.names=colnames(x.variates)
    temp.samps=rownames(x.variates)
    x.variates=as.data.frame(cbind(-1*x.variates[,1],x.variates[,-1]))
    colnames(x.variates)=temp.names
    rownames(x.variates)=temp.samps
    
    temp.names = colnames(x.loadings)
    temp.samps = rownames(x.loadings)
    x.loadings = as.data.frame(cbind(-1 * x.loadings[, 1], 
                                     x.loadings[, -1]))
    colnames(x.loadings) = temp.names
    rownames(x.loadings) = temp.samps
    
    
  }

  
  name=sub(".txt","",file)
  savename=paste(name,"_prcomp_scores.txt",sep='');
  write.table(cbind(Score=rownames(x.variates),x.variates),savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_loadings.txt",sep='');
  write.table(cbind(Loadings=rownames(x.loadings),x.loadings),savename,sep='\t',row.names=FALSE,quote=FALSE);
  savename=paste(name,"_prcomp_sdev.txt",sep='');
  write.table(pca_evalues,savename,sep='\t',row.names=FALSE,quote=FALSE);
  print(summary(pca))
  screeplot(pca)
  
 
  
  
  x.variates$type = sample.type[match(rownames(x.variates), 
                                      sample.names)]
  
  
  
  
  pc.pred = ggplot(data = x.variates, aes_string(x = comp.x, 
                                                 y = comp.y)) + geom_point(size = I(2), aes(color = factor(type))) + 
    theme(legend.position = "right", plot.title = element_text(size = 30), 
          legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
          axis.title = element_text(size = 30), legend.background = element_rect(), 
          axis.text.x = element_text(margin = margin(b = -2)), 
          axis.text.y = element_text(margin = margin(l = -14))) + 
    labs(title = title) + theme_bw() + if (labels == TRUE) {
      geom_text(data = x.variates, mapping = aes(label = (rownames(x.variates))), 
                check_overlap = TRUE, size = 2.5)
    }
  pc.pred
  ggsave(paste0(output_folder, train_string, ".PCA.",comp.x, "_vs_", comp.y, savetype), 
         dpi = 300, plot = pc.pred, width = w, height = h)
  
  ### KEEP IT THIS WAY. 
  if (varimax == T) {
    temp_names=colnames(rotated.data2)
    temp.samps=rownames(rotated.data2)
    rotated.data2 <- as.matrix(rotated.data2[, 1:(varimax.comp)]) %*%
      rotation$rotmat
    colnames(rotated.data2) = temp_names[1:(varimax.comp)]
    rownames(rotated.data2) = temp.samps
    rotated.data2 = as.data.frame(rotated.data2)
  }

  if(rotate==T){
    temp.names=colnames(rotated.data2)
    temp.samps=rownames(rotated.data2)
    rotated.data2=as.data.frame(cbind(-1*rotated.data2[,1],rotated.data2[,-1]))
    colnames(rotated.data2)=temp.names
    rownames(rotated.data2)=temp.samps
  }
  
  rotated.data2=as.data.frame(rotated.data2)
  write.table(cbind(Sample = rownames(rotated.data2), (rotated.data2)), 
              paste0(output_folder, test_string, "_projected_onto_", 
                     train_string, "_PCA_predicted.scores.txt"), sep = "\t", 
              row.names = F, quote = F)
  
  prediction.to.write=scale(rotated.data2)
  prediction.to.write=as.data.frame(prediction.to.write)
  prediction.to.write$type = sample.type[match(rownames(prediction.to.write),sample.names)]
  prediction.to.write$prediction =ifelse(prediction.to.write$PC1 >=threshold, 1,0)
  
  write.table(prediction.to.write,paste0(output_folder,test_string,"_projected_onto_",train_string,"_",comps,"_comps_PCA_prediction.txt"),col.names=NA,quote=F,sep="\t",row.names=T)
  
  rotated.data2$type = sample.type[match(rownames(rotated.data2),sample.names)]
  
  
  
  pc.pred2 <- ggplot(rotated.data2, aes_string(x = comp.x, y = comp.y)) + 
    geom_point(size = I(2), aes(color = factor(type)))
  
  if(do.legend==F) {
    
    
    pc.pred2<- pc.pred2 +geom_point(size = I(2), show.legend=F,aes(color = factor(type)))+ theme_bw()  + theme( legend.position = "none",plot.title = element_text(size = 30), 
                               axis.title = element_text(size = 30), legend.background = element_rect(), 
                               axis.text.x = element_text(margin = margin(b = -2)), 
                               axis.text.y = element_text(margin = margin(l = -14))) + 
      labs(title = title) + guides(fill=FALSE)+
     if (labels == TRUE) {
        geom_text(data = prediction, mapping = aes(label = (rownames(prediction))), 
                  check_overlap = TRUE, size = 2.3)
      }} else{
        pc.pred2<- pc.pred2 +theme(legend.position = "right", plot.title = element_text(size = 30), 
                                   legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
                                   axis.title = element_text(size = 30), legend.background = element_rect(), 
                                   axis.text.x = element_text(margin = margin(b = -2)), 
                                   axis.text.y = element_text(margin = margin(l = -14))) + 
          guides(color = guide_legend(title = "Type")) + labs(title = title) + 
          theme_bw()  + if (labels == TRUE) {
            geom_text(data = prediction, mapping = aes(label = (rownames(prediction))), 
                      check_overlap = TRUE, size = 2.3)
          }}
  pc.pred2
  
  
  if (plot_both == T) {
    comb = rbind(rotated.data2, x.variates)
    pc.pred3 = ggplot(data = comb, aes_string(x = comp.x, 
                                              y = comp.y), ) + geom_point(size = I(3), aes(color = factor(type), 
                                                                                           shape = factor(type))) + theme(legend.position = "right", 
                                                                                                                          plot.title = element_text(size = 30), legend.text = element_text(size = 22), 
                                                                                                                          legend.title = element_text(size = 20), axis.title = element_text(size = 30), 
                                                                                                                          legend.background = element_rect(), axis.text.x = element_text(margin = margin(b = -2)), 
                                                                                                                          axis.text.y = element_text(margin = margin(l = -14))) + labs(title = title) + theme_bw() + 
      if (labels ==  TRUE) {
        geom_text(data = comb, mapping = aes(label = (rownames(comb))), 
                  check_overlap = TRUE, size = 2.5)
      }
    if (!is.null(shape.palette)) {
      pc.pred3 <- pc.pred3 + scale_shape_manual(legendname, 
                                                values = shape.palette)
    }
    if (!is.null(colpalette)) {
      pc.pred3 <- pc.pred3 + scale_color_manual(legendname, 
                                                values = colpalette)
    }
    if (ellipses == T) {
      pc.pred3 <- pc.pred3 + stat_ellipse(aes(color = factor(type)), 
                                          level = conf)
    }
    if(do.legend==F){
      pc.pred3<-pc.pred3 +theme(legend.position="none")
    }
    if (saveplot == T) {
      ggsave(paste0(output_folder, test_string, "_projected_onto_", 
                    train_string, "_", comp.x, "_vs_", comp.y, savetype), 
             dpi = 300, plot = pc.pred3, width = w, height = h)
    }
    pc.pred3
  }
  else {
    pc.pred = pc.pred + geom_point(data = x.variates, aes_string(x = comp.x, 
                                                                 y = comp.y)) + geom_point(size = I(1.3), aes(color = factor(type))) + 
      theme(legend.position = "right", plot.title = element_text(size = 30), 
            legend.text = element_text(size = 22), legend.title = element_text(size = 20), 
            axis.title = element_text(size = 30), legend.background = element_rect(), 
            axis.text.x = element_text(margin = margin(b = -2)), 
            axis.text.y = element_text(margin = margin(l = -14))) + 
      labs(title = title) + theme_bw()
    if (saveplot == T) {
      ggsave(paste0(output_folder, test_string, "_projected_onto_", 
                    train_string, "_", comp.x, "_vs_", comp.y, savetype), 
             dpi = 300, plot = pc.pred, width = w, height = h)
    }
    pc.pred
  }
