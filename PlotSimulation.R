library(ggplot2)

Plot <- function(result,type="all"){
  
  df=data.frame(as.numeric(result[,1]))
  df=cbind.data.frame(df,as.numeric(result[,2]))
  df=cbind.data.frame(df,result[,3])
  df=cbind.data.frame(df,c('0.MAR','3.Model','3.Model','3.Model','3.Model','0.MAR','0.MAR','2.Implicit','2.Implicit','0.MAR','0.MAR','2.Implicit','2.Implicit','2.Implicit','2.Implicit','0.MAR','0.MAR','2.Implicit','2.Implicit','0.MAR','2.Implicit'))
  df=cbind.data.frame(df,c('Mean','Pred','Pred','Tot','Tot','Tot','Pred','Tot','Pred','Tot','Pred','Tot','Pred','Tot','Pred','Tot','Pred','Tot','Pred','Mean','Mean'))
  colnames(df)=c("mse","mseTrue","algo","meth","type") 
  df1=df[df$type=='Tot' | df$type=='Mean',]
  df1$algo=rep(c('1.mean','3.soft','4.FISTA','3.soft','3.soft','2.PCA','2.PCA','5.mimi','4.FISTA','4.FISTA','6.rf','6.rf'),50)
  df1=df1[df1$algo!='6.rf',]
  if (type=="WithoutMask"){
    df1=df1[df1$meth!='2.Implicit',]
  }
  
  df2=df[df$type=='Pred' | df$type=='Mean',]
  df2$algo=rep(c('1.mean','3.soft','4.FISTA','3.soft','3.soft','2.PCA','2.PCA','5.mimi','4.FISTA','4.FISTA','6.rf','6.rf'),50)
  df2=df2[df2$algo!='6.rf',]
  if (type=="WithoutMask"){
    df2=df2[df2$meth!='2.Implicit',]
  }
  
  
  get_legend<-function(myggplot){
    tmp <- ggplot_gtable(ggplot_build(myggplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  if (type=="all"){
    vec_color <- c("#0033FF","#006600","#CC0000")
  }else if (type=="WithoutMask"){
    vec_color <- c("#0033FF","#CC0000")
  }
  plot1<-ggplot(df1)+aes(x=algo,y=mseTrue,colour=meth)+geom_boxplot()+xlab("")+ylab("")+scale_colour_manual('',labels=c('MAR','Mask','Model'),values=vec_color)+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),axis.text.x=element_text(size=24, angle=45),axis.text.y=element_text(size=24),title=element_text(size=24),legend.text=element_text(size=24))+ labs(fill = "")+scale_x_discrete(label=c("mean","PCA","soft","FISTA","mimi")) #+scale_y_continuous(limits=c(0.5,1.6))
  legend <- get_legend(plot1)
  plot1<-ggplot(df1)+aes(x=algo,y=mseTrue,colour=meth)+geom_boxplot()+xlab("")+ylab("")+scale_colour_manual('meth',values=vec_color)+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),axis.text.x=element_text(size=24, angle=45),axis.text.y=element_text(size=24),title=element_text(size=24),legend.text=element_text(size=24))+theme(legend.position="none")+scale_x_discrete(label=c("mean","PCA","soft","FISTA","mimi")) #+scale_y_continuous(limits=c(0.5,1.6))
  
  plot2<-ggplot(df2)+aes(x=algo,y=mse,colour=meth)+geom_boxplot()+xlab("")+ylab("")+scale_colour_manual('meth',values=vec_color)+theme(axis.title.x=element_text(size=24),axis.title.y=element_text(size=24),axis.text.x=element_text(size=24, angle=45),axis.text.y=element_text(size=24),title=element_text(size=24),legend.text=element_text(size=24))+theme(legend.position="none")+scale_x_discrete(label=c("mean","PCA","soft","FISTA","mimi")) #+scale_y_continuous(limits=c(0.5,1.6))
  
  library(gridExtra)
  grid.arrange(plot1,plot2,legend, ncol=3, widths=c(2.4, 2.4, 0.6))

}


