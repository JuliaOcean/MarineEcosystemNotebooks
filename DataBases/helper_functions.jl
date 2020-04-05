
module cmap_helpers

using Pandas, PyCall
PyCmap = pyimport("pycmap")
#cmap = PyCmap.API(token="your-own-API-key")
cmap = PyCmap.API()

"""
    get(t::String,v::String)

Retrieve variable v from CMAP table t. Return along with
meta-data, position and time, in a Dict.
"""
function get(t::String,v::String)
    df=Pandas.DataFrame(cmap.get_dataset(t))
    me=Pandas.DataFrame(cmap.get_metadata(t,v))
    x=Dict("Variable" => v,"Unit" => deepcopy(values(me["Unit"])[1]),
        "Long_Name" => deepcopy(values(me["Long_Name"])[1]),
        "Data_Source" => deepcopy(values(me["Data_Source"])[1]),
        "lon" => deepcopy(values(df[:lon])),
        "lat" => deepcopy(values(df[:lat])),
        "time" => deepcopy(values(df[:time])),
        "val" => deepcopy(values(df[Symbol(v)])))
    return x
end

"""
    tables(ListName::String)

List of CMAP tables for subset `ListName`.
"""
function tables(ListName::String)
    list0=[];
    if ListName=="G3"
        list0=["tblKM1906_Gradients3","tblKM1906_Gradients3_uway_optics","tblKM1906_Gradients3_uwayCTD","tblKM1906_Gradients3_uw_tsg"]
    elseif ListName=="G2"
        list0=["tblMGL1704_Gradients2_CTD","tblMGL1704_Gradients2_uway_optics"]
    elseif ListName=="G1"
        list0=["tblKOK1606_Gradients1_CTD","tblKOK1606_Gradients1_uway_optics"]
    elseif ListName=="MoreG1"
        list0=["tblKOK1606_Gradients1_Nutrients","tblKOK1606_Gradients1_Dissolved_Gasses",
        "tblKOK1606_Gradients1_TargetedMetabolites","tblKOK1606_Gradients1_Diazotroph"]
    elseif ListListNameId=="MoreG2"
        list0=["tblMGL1704_Gradients2_Nutrients","tblMGL1704_Gradients2_Diazotroph",
        "tblMGL1704_Gradients2_TargetedMetabolites","tblMGL1704_Gradients2_Trace_Metals"]
    end
    return list0
end

end
