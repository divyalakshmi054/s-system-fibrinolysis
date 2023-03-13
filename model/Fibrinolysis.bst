// S-law model of fibrinolysis
// Author: dl653@cornell.edu

#static::start
    tPA,PAI1,TAFI,AT
#static::end

#dynamic::start
    FII,FI,FIIa,FIa,FDP,Plgn,Plasmin,FIIa_inactive
#dynamic::end

#structure::start
    // name::{reactants} --> {products}
    // 1. name: reaction name (used in rate section)
    // 2. {reactants}: list of dynamic species acting as a reactants
    // 3. {products}: list of dynamic species acting as product
    
    // reactions in the pathway -
    r1::{FII} --> {FIIa}
    r2::{FIIa} --> {FIIa_inactive}
    r3::{FI} --> {FIa}
    r4::{Plgn} --> {Plasmin}
    r5::{FIa} --> {FDP}
#structure::end

#rate::start
    // name::{factor_1,factor_2,....}
    // 1: name = reaction name (must be in the #pragma::structure section)
    // 2: list of factor symbols (can be both static and dynamic, order does not matter)
    // 3: the order of the records does not matter here, we establish the order from the structure section
    r1::{FII,FIIa}
    r2::{FIIa,AT}
    r3::{FI,FIIa}
    r4::{Plgn,tPA,PAI1}
    r5::{FIa,Plasmin,TAFI}
#rate::end