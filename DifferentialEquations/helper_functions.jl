
"""
    SetupPeriodicDomain(np::Integer=16)

Set up a periodic domain of size np x np

```
np=16 #domain size is np x np
Γ=SetPeriodicDomain(np)
```
"""
function SetupPeriodicDomain(np::Integer=16)
    γ,Γ=GridOfOnes("PeriodicDomain",1,np)
    Γ["XC"][1]=vec(0.5:1.0:np-0.5)*ones(1,np)
    Γ["XG"][1]=vec(0.0:1.0:np-1.0)*ones(1,np)
    Γ["YC"][1]=ones(np,1)*transpose(vec(0.5:1.0:np-0.5))
    Γ["YG"][1]=ones(np,1)*transpose(vec(0.0:1.0:np-1.0))
    return Γ
end

