
using MeshArrays, Plots, CSV, DataFrames, JLD

GCMGridSpec("LLC90")
GCMGridLoad()

"""
    CreateModelSamples(iFile,iVar,indx)

Example:
```
df, mp = CreateModelSamples(2,1,[])

df2, mp2 = CreateModelSamples(2,2,df.indx)
df = join(df,df2,on = [:indx,:lon,:lat])
mp = merge(mp,mp2)

save("mp.jld", mp)
mp2=load("mp.jld")

CSV.write("df.csv", df)
df2 = CSV.File("df.csv") |> DataFrame!
```
"""
function CreateModelSamples(iFile,iVar,indx)

    nSample=100

    #test cases
    if (iFile==1)
        InDir="occci-daily-llc90/"
        InFil="OC_CCI_L3S_Rrs_412_2010"
        InPrec=Float32
        varName="Rrs412"
    else #minimal test case
        InDir="GRID_LLC90/"
        InFil="Depth.data"
        InPrec=Float64
        varName="Depth"
    end
    #... later this should also treat model output
    #... and variable lists -> 1 csv file + n binaries

    #here should be a loops over
    #1) days
    #2) variables

    tmp=read_bin(InDir*InFil,InPrec)
    m=fill(1., MeshArrays.XC) + 0. * mask(MeshArrays.hFacC[:,:,1],NaN,0)
    tmp=m*tmp[:,:,1]
    tmp=mask(tmp,NaN,0)

    #extract random subset for each day
    lon=convert2gcmfaces(MeshArrays.XC)
    lat=convert2gcmfaces(MeshArrays.YC)
    df=convert2gcmfaces(tmp)
    isempty(indx) ? indx=rand(findall(df.>0.),nSample) : nothing
    #... store in a DataFrame(Row)
    #df = DataFrame(indx = indx, lon = lon[indx], lat = lat[indx], x = tmp[indx])
    df = DataFrame( [indx, lon[indx], lat[indx], df[indx] ],[:indx, :lon, :lat, Symbol(varName)])

    #get full fields for january 1
    mp=convert2array(tmp)
    mp=mp[:,62:239]
    mp=Dict( varName => circshift(mp,(142,0)))
    #... store to file
    #heatmap(tmp)

    #... add dataframe(row) to dataframe
    #... close loops

    #... write dataframe to file
    return df, mp

end
