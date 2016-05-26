library(CRHMr)

obs<-readObsFile('Marmot_Hourly_ArrayMetData_withT_g_1Oct05-30Sept14_update_19Nov2014.obs',timezone = 'CST')

#Marmot Creek basin array data from
# alpine = 1
#fisera ridge =2 
#vista view = 3
#upper clearing = 4
#upper forest =5
#hay meadow =6
#level forest staions =7

# t 7 (°C)																																							
# rh 7 (%)
# T_g 7 (°C)
# u 7 (m/s)																																							
# p 7 (mm/int)																																							
# Qsi 7 (W/m^2)

fr<-obs[c('datetime','t.2','rh.2','u.2','p.2','Qsi.2')]
fr$datetime<-format(fr$datetime,format = '%Y%m%dT%H%M%S')
names(fr)<-c('datetime','t','rh','u','p','Qsi')
write.table(fr,file="fr_2005_2014.txt",sep="\t",quote=F,row.names=F)


vv<-obs[c('datetime','t.3','rh.3','u.3','p.3','Qsi.3')]
vv$datetime<-format(vv$datetime,format = '%Y%m%dT%H%M%S')
names(vv)<-c('datetime','t','rh','u','p','Qsi')
write.table(vv,file="vv_2005_2014.txt",sep="\t",quote=F,row.names=F)

uc<-obs[c('datetime','t.4','rh.4','u.4','p.4','Qsi.4')]
uc$datetime<-format(uc$datetime,format = '%Y%m%dT%H%M%S')
names(uc)<-c('datetime','t','rh','u','p','Qsi')
write.table(uc,file="uc_2005_2014.txt",sep="\t",quote=F,row.names=F)

hm<-obs[c('datetime','t.6','rh.6','u.6','p.6','Qsi.6')]
hm$datetime<-format(hm$datetime,format = '%Y%m%dT%H%M%S')
names(hm)<-c('datetime','t','rh','u','p','Qsi')
write.table(hm,file="hm_2005_2014.txt",sep="\t",quote=F,row.names=F)
