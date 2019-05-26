
using MeshArrays, Plots

GCMGridSpec("LLC90")
GCMGridLoad()

#test cases
if false
    InDir="occci-daily-llc90/"
    InFil="OC_CCI_L3S_Rrs_412_2010"
    InPrec=Float32
else #minimal test case
    InDir="GRID_LLC90/"
    InFil="Depth.data"
    InPrec=Float64
end
#... later this should also treat model output
#... and variable lists -> 1 csv file + n binaries

#here should be a loops over
#1) days
#2) variables

D=read_bin(InDir*InFil,InPrec)

#extract random subset for each day
tmp=convert2gcmfaces(D[:,:,1])
ii=findall(tmp.>0.)
jj=rand(ii,100)
tmp[jj]
#... store in a DataFrame(Row)

#get full fields for january 1
#tmp=convert2array(MeshArrays.XC)
tmp=convert2array(D[:,:,1])
tmp=tmp[:,62:239]
tmp=circshift(tmp,(142,0))
#... store to file
heatmap(tmp)

#... add dataframe(row) to dataframe
#... close loops

#... write dataframe to file
