library(CRHMr)


obs<-readObsFile('Marmot_Hourly_ArrayMetData_withT_g_1Oct05-30Sept17_update_10Jan2018.obs',timezone = 'etc/GMT+7')


#Marmot Creek basin array data from
# alpine = 1
#fisera ridge =2 
#vista view = 3
#upper clearing = 4
#upper forest =5
#hay meadow =6
#level forest staions =7

# fr<-obs[c('datetime','t.2','rh.2','u.2','p.2','Qsi.2')]
# fr$datetime<-format(fr$datetime,format = '%Y%m%dT%H%M%S')
# names(fr)<-c('datetime','t','rh','u','p','Qsi')
# write.table(fr,file="fr_2005_2014.txt",sep="\t",quote=F,row.names=F)
# 
# 
# vv<-obs[c('datetime','t.3','rh.3','u.3','p.3','Qsi.3')]
# vv$datetime<-format(vv$datetime,format = '%Y%m%dT%H%M%S')
# names(vv)<-c('datetime','t','rh','u','p','Qsi')
# write.table(vv,file="vv_2005_2014.txt",sep="\t",quote=F,row.names=F)

uc<-obs[c('datetime','t.4','rh.4','u.4','p.4','Qsi.4','T_g.4')]
names(uc)<-c('datetime','t','rh','u','p','Qsi','T_g')
uc$hydro <- hydroYear(uc)
uc[uc['u']>4 & uc['hydro']==2013,'u'] <- 4
uc$hydro <- NULL #drop

uc$datetime<-format(uc$datetime,format = '%Y%m%dT%H%M%S')
write.table(uc,file="uc_2005_2018.txt",sep="\t",quote=F,row.names=F)

# hm<-obs[c('datetime','t.6','rh.6','u.6','p.6','Qsi.6')]
# hm$datetime<-format(hm$datetime,format = '%Y%m%dT%H%M%S')
# names(hm)<-c('datetime','t','rh','u','p','Qsi')
# write.table(hm,file="hm_2005_2014.txt",sep="\t",quote=F,row.names=F)


    
                   
                   