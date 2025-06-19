# Importing and pre-processing spectra, and predicting foliar traits using 
# pre-trained PLSR models

# loading libraries
library(pacman)
p_load(spectrolab, doParallel, tidyr, cowplot, reshape2, stringr, tools, GGally)

# Importing spectra
set.seed(112345)
cores=detectCores()-1

split_path=function(path) {
  if (dirname(path) %in% c("_", path)) return(basename(path))
  return(c(basename(path),split_path(dirname(path))))
}

cn=function(x)
{
  print(as.matrix(colnames(x)))
}

for (j in c("asd","sed","sig"))
{
  if (exists('finaloutput')){rm(finaloutput)}
  folderlist=unique(sub("/[^/]+$","",list.files(path="Demo/Spectra_test",
                                                pattern=paste0(".",j,"$"),
                                                recursive=TRUE,
                                                full.names = TRUE)))
  if (length(folderlist)!=0){
    cl=makeCluster(cores)
    registerDoParallel(cl)
    
    finaloutput=foreach(i=folderlist,
                        .combine=rbind,
                        .packages=c('spectrolab'))%dopar%
      {
        spectra=read_spectra(path=i,format=j)
        data.frame(spectra,folder=i,format=j,check.names=FALSE)
      };FST=as_spectra(finaloutput[,1:2152],name_idx=1)
    
    stopCluster(cl)
    
    if (j=="asd"){
      FST=match_sensors(
        FST,
        splice_at=c(1001,1801)
      )
    }
    if (j=="sig"){
      FST=match_sensors(
        FST,
        splice_at=c(990,1900)
      )
    };FSTdf=data.frame(FST,check.names=FALSE)
    if (!exists('FFO')){FFO=FSTdf[0,]}
    FFO=rbind(FFO,FSTdf)
  }};FSTfinal=as_spectra(FFO,name_idx=1);FSTfinal_norm=normalize(FSTfinal);FSTfinal_sd=sd(FSTfinal)

FST_f=data.frame(FSTfinal)
FST_n=data.frame(FSTfinal_norm)

head(FST_f)[c(1,402:412)] #Original Spectra

#----Process Spectra, Vector Normalize, and Resample----
headcount=ncol(FST_f)-2151
specstart=headcount+1

cl=makeCluster(cores);registerDoParallel(cl);sampledata_VN=foreach(i=1:nrow(FST_f),
                                                                   .combine=rbind)%dopar%
  {
    a=FST_f[i,]
    b=sqrt(rowSums(a[,specstart:ncol(a)]^2))
    d=a[,-c(1:headcount)]
    sample_name=a[,c(1:headcount)]
    f=d/b
    cbind(sample_name,f)
  }
stopCluster(cl)

head(sampledata_VN)[c(1,402:412)] #Vector Normalized Spectra

sampledata=sampledata_VN

sampledata_wav=sampledata[,-c(1:headcount)]
sampledata_head=sampledata[,c(1:headcount)]
sampledata_5nm_wav=sampledata_wav[,seq_len(ncol(sampledata_wav)) %% 5 == 1]
sampledata_5nm=cbind(sampledata_head,sampledata_5nm_wav)
sampledata=sampledata_5nm

head(sampledata)[c(1,82:92)] #Vector Normalized Resampled Spectra

# Plot spectra

plot(FSTfinal_norm, ylab = "Vector Normalized Reflectance %", xlab = "Wavelength (nm)", lty = 1, type = "l")

plot_regions(
  FSTfinal,
  regions = default_spec_regions(),
  col = grDevices::rgb(0.7, 0.7, 0.7, 0.3),
  border = FALSE,
  add = TRUE,
  add_label = TRUE,
  cex_label = 1,
)

#Outlier Detection
minmin=which.min(FST_f$X1200)
min(FST_f$X1200,na.rm=T)
FST_f[minmin,1]

#Outlier Detection
maxmax=which.max(FST_f$X550)
max(FST_f$X550,na.rm=T)
FST_f[maxmax,1]

# #Correlation Plot
ggcorr(sampledata[c(-1)])

# Grouped Outlier Detection Method
sample_data_plot=FST_f

sample_data_plot$Sample=str_split_fixed(sample_data_plot$sample_name,"_",n=4)[,1]
sample_data_plot$Replicate=str_split_fixed(sample_data_plot$sample_name,"_",n=4)[,2]

sample_data_plot_sub=sample_data_plot[,c(ncol(sample_data_plot)-1,ncol(sample_data_plot),4:ncol(sample_data_plot)-2)]

tic=0;for (i in unique(sample_data_plot_sub$Sample)){
  a=subset(sample_data_plot_sub,sample_data_plot_sub$Sample==i)
  aa=a[,-c(1:2)]
  b=sqrt(apply(aa[,1:2151]^2,1,sum))
  d=aa[,1:2151]/b
  e=cbind(a[,1:2],d)
  f=melt(e, id.vars=c("Sample","Replicate"))
  f$Replicate=as.factor(f$Replicate)
  print(ggplot(f, aes(x=as.numeric(variable), y=value, group=Replicate,color=Replicate))+
          geom_line()+
          ggtitle(paste("Sample:",f$Sample[1]))+
          ylab("value")+
          theme_bw()+
          scale_x_continuous(name="band",
                             limits=c(0, 2151),
                             breaks=seq(150, 2151,by=500),
                             labels=seq(500,2500, by=500)))
  tic=tic+1;print(paste0(tic,"/",length(unique(sample_data_plot_sub$Sample))))
  readline(prompt="Press [enter] to continue")
};rm(a,b,d,e,f,i,sample_data_plot_sub,tic)

#----Predict Traits----

modeldirectory=list.files(path="Models/models_PSR_fresh_spectra",
                          pattern=".csv$",
                          recursive=FALSE,
                          full.names = TRUE)

cl=makeCluster(cores);registerDoParallel(cl);finished_output=foreach(model=modeldirectory,
                                                                     .combine=rbind)%do%
  {
    modelname=split_path(model)[1]
    modelcoeffs_mat=data.matrix(read.csv(model, header = TRUE, row.names = 1,sep = ","))
    matplottable=modelcoeffs_mat[1,c(12:412)]
    matplot(matplottable, main=modelname, type = 'l', xaxt="n", xlab="Wavelength (nm)", ylab="Coefficient")
    axis(1, at=c(seq(from=0,to=400,by=2)), labels=c(seq(from=400,to=2400,by=10)))
    sampledata_tail=sampledata[,-c(1:headcount)]
    sampledata_head=sampledata[,c(1:headcount)]
    sampledata_tail$constant=1
    sampledata_tail=cbind(intercept=sampledata_tail[,ncol(sampledata_tail)],sampledata_tail[,c(1:ncol(sampledata_tail)-1)])
    sampledata_tail_mat=as.matrix(sampledata_tail)
    modeloutput=sampledata_tail_mat %*% t(modelcoeffs_mat)
    modeloutput_bind=cbind(sampledata_head,as.data.frame(modeloutput))
    data.frame(modelname=modelname,
               sample_name=sampledata_head,
               t_mean=rowMeans(modeloutput),
               t_std=apply(modeloutput,1,sd),
               as.data.frame(modeloutput),
               check.names=FALSE)
  };stopCluster(cl)

predplot=finished_output[,c(1,3:4)]
predplot$Group=str_split_fixed(finished_output$sample_name,"_",n=2)[,1]

# Add species names as a column in the dataframe for plotting
custom_labs <- c("Brun" = "Berzelia lanuginosa",
                 "ELFI" = "Elegia cuspidata",
                 "LECO" = "Leecadendron coniferum",
                 "T068" = "Elegia filacea",
                 "T081" = "L. lau (yellow leaves)",
                 "T085" = "Erica spp",
                 "T68T" = "E. fil (inflorescence)",
                 "T81N" = "Leucadendron laureolum",
                 "T85T" = "Erica spp (flowers)")
predplot$custom_labs <- custom_labs[predplot$Group]

predsub=finished_output[,c(1:4)]

predspread=pivot_wider(data=predsub,id_cols=sample_name,names_from=modelname,values_from=c("t_mean","t_std"))
predunlist=as.data.frame(unnest(predspread))
colnames(predunlist)=sub(".csv", "", colnames(predunlist))
# write.csv(predunlist,"SpectroPredictR_finished_model_output.csv",row.names=FALSE)
predunlist[1:4,]

# Plot Boxplots 
for (i in unique(predplot$modelname)){
  
  datanow=subset(predplot,modelname==i)
  
  d=ggplot(datanow,
           aes(x=custom_labs,y=t_mean))+
    ggtitle(paste0(file_path_sans_ext(i)," (Fresh Samples) [n=",nrow(datanow),"]"))+
    geom_boxplot(color="blue",fill="lightblue",size=0.5)+
    coord_flip()+
    xlab("Sample")+
    ylab("Mean of Predictions")+
    scale_color_brewer(palette="Set1")+
    theme(text = element_text(size=15))+
    scale_x_discrete(limits=rev(levels(as.factor(datanow$custom_labs))))
  
  e=ggplot(datanow,
           aes(x=custom_labs,y=t_std,))+
    ggtitle("    ")+
    geom_boxplot(color="black",fill="yellow",size=0.5)+
    ylab("StDev of Predictions")+
    scale_color_brewer(palette="Set1")+
    coord_flip()+
    scale_x_discrete(limits=rev(levels(as.factor(datanow$custom_labs))))+
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none")
  
  print(plot_grid(d,e,rel_widths = c(4,3)))
  
  ggsave(
    file=paste0("visuals/",file_path_sans_ext(i),"_fresh_plot.png"),
    width = 15,
    height = 8.5,
    units="in",
    dpi = 300
  )
}
