with(LinearAlgebra):

JacX := Matrix(42, 42, [[0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[4631863,0,5210176,-771100,-12225832,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0],[4541169434870,0,0,23264218531158,24069354941916,0,0,0,0,0,0,0,-771100,0,0,0,0,0,0,0,-12225832,0,0,0,0,0,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0],[0,-6174063,4644266,5402963,-2371390,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-771100,3973466,-5402963,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,739882,0,0,485360,7015656],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0],[-367266045257292086464,0,0,-161832119336657612096,66540238208918325394,0,0,0,0,0,0,-771100,46528437062316,0,0,0,0,0,0,-12225832,48138709883832,0,0,0,0,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,41987267627446,-24069354941916,-18723049096288,0,0,0,0,0,0,0,0,5402963,0,0,0,0,0,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,23264218531158,0,18723049096288,0,0,0,0,0,0,1,3973466,0,0,0,0,0,0,0,-5402963,0,0,0,0,0,0,0,0,-771100,0,0,0,0,0,0,739882,0,0,0,29796571083828,-24069354941916],[0,0,-771100,-4243028,-5402963,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,8524144,485360,-7015656,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[3321411619153482393725922562,0,0,-796067883340498923546834134,-2218679076696712477258179848,0,0,0,0,0,-771100,69792655593474,-485496358009972836288,0,0,0,0,0,-12225832,72208064825748,199620714626754976182,0,0,0,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,43601806583976862272,-66540238208918325394,-205433925920634474368,0,0,0,0,0,0,0,5402963,-37446098192576,0,0,0,0,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,-48138709883832,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-161832119336657612096,0,205433925920634474368,0,0,0,0,0,1,3973466,0,0,0,0,0,0,0,-5402963,37446098192576,0,0,0,0,0,0,0,-771100,46528437062316,0,0,0,0,0,739882,0,0,0,0,-496250264691182602834,-66540238208918325394],[0,0,23264218531158,0,18723049096288,0,0,0,0,0,0,0,-4243028,0,0,0,0,0,0,0,-5402963,0,0,0,0,0,0,0,0,-771100,0,0,0,0,0,1,8524144,0,29796571083828,24069354941916,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[30742876800843791320167620017042920,0,0,47793818282399302637313018375608180,187093231355119609259230575406078,0,0,0,0,-771100,93056874124632,-970992716019945672576,-3184271533361995694187336536,0,0,0,0,-12225832,96277419767664,399241429253509952364,-8874716306786849909032719392,0,0,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-4913547385834480240819590830,2218679076696712477258179848,4117479502493981317272756696,0,0,0,0,0,0,5402963,-56169147288864,-616301777761903423104,0,0,0,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,-72208064825748,-199620714626754976182,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-796067883340498923546834134,0,-4117479502493981317272756696,0,0,0,0,1,3973466,0,0,0,0,0,0,0,-5402963,56169147288864,616301777761903423104,0,0,0,0,0,0,-771100,69792655593474,-485496358009972836288,0,0,0,0,739882,0,0,0,0,0,5534893686620078780470129720,2218679076696712477258179848],[0,0,-161832119336657612096,0,205433925920634474368,0,0,0,0,0,0,-4243028,0,0,0,0,0,0,0,-5402963,37446098192576,0,0,0,0,0,0,0,-771100,46528437062316,0,0,0,0,1,8524144,0,0,-496250264691182602834,66540238208918325394,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-1290768342652267516443479931546023286146818,0,0,-447426500232890640457842369370436667327162,894358278459026413928193823774198568360308,0,0,0,-771100,116321092655790,-1618321193366576120960,-7960678833404989235468341340,238969091411996513186565091878040900,0,0,0,-12225832,120346774709580,665402382089183253940,-22186790766967124772581798480,935466156775598046296152877030390,0,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,64844759763954813954458416734173440,-187093231355119609259230575406078,-17050941481555511317145398358565260,0,0,0,0,0,5402963,-74892196385152,-1232603555523806846208,16469918009975925269091026784,0,0,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,-96277419767664,-399241429253509952364,8874716306786849909032719392,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,47793818282399302637313018375608180,0,17050941481555511317145398358565260,0,0,0,1,3973466,0,0,0,0,0,0,0,-5402963,74892196385152,1232603555523806846208,-16469918009975925269091026784,0,0,0,0,0,-771100,93056874124632,-970992716019945672576,-3184271533361995694187336536,0,0,0,739882,0,0,0,0,0,0,-25042392508547629407962232317859806,-187093231355119609259230575406078],[0,0,-796067883340498923546834134,0,-4117479502493981317272756696,0,0,0,0,0,-4243028,0,0,0,0,0,0,0,-5402963,56169147288864,616301777761903423104,0,0,0,0,0,0,-771100,69792655593474,-485496358009972836288,0,0,0,1,8524144,0,0,0,5534893686620078780470129720,-2218679076696712477258179848,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[2970834789581641652400857615710374522053581511540,0,0,-15038070868923163890624777882201038403033277340040,-17325082089263041613534300922518663065279234311750,0,0,-771100,139585311186948,-2427481790049864181440,-15921357666809978470936682680,716907274235989539559695275634122700,-2684559001397343842747054216222620003962972,0,0,-12225832,144416129651496,998103573133774880910,-44373581533934249545163596960,2806398470326794138888458631091170,5366149670754158483569162942645191410161848,0,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,395915342186486235527795192805149951492494,-894358278459026413928193823774198568360308,-843341842419376875985637562175586618819656,0,0,0,0,5402963,-93615245481440,-2054339259206344743680,41174795024939813172727566960,-85254707407777556585726991792826300,0,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,-120346774709580,-665402382089183253940,22186790766967124772581798480,-935466156775598046296152877030390,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-447426500232890640457842369370436667327162,0,843341842419376875985637562175586618819656,0,0,1,3973466,0,0,0,0,0,0,0,-5402963,93615245481440,2054339259206344743680,-41174795024939813172727566960,85254707407777556585726991792826300,0,0,0,0,-771100,116321092655790,-1618321193366576120960,-7960678833404989235468341340,238969091411996513186565091878040900,0,0,739882,0,0,0,0,0,0,0,-700958984478524379276471733774551146637668,-894358278459026413928193823774198568360308],[0,0,47793818282399302637313018375608180,0,17050941481555511317145398358565260,0,0,0,0,-4243028,0,0,0,0,0,0,0,-5402963,74892196385152,1232603555523806846208,-16469918009975925269091026784,0,0,0,0,0,-771100,93056874124632,-970992716019945672576,-3184271533361995694187336536,0,0,1,8524144,0,0,0,0,-25042392508547629407962232317859806,187093231355119609259230575406078,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[810328895419757290411394963849398487057516046341935956802,0,0,564149090960036788255909059719299671527821927753365486850,-253774647033375172525007962942617633458863270496080865904,0,-771100,162849529718106,-3398474506069809854016,-27862375916917462324139194690,1672783639883975592305955643146286300,-9395956504890703449614689756779170013870402,-105266496082462147234373445175407268821232941380280,0,-12225832,168485484593412,1397345002387284833274,-77653767684384936704036294680,6548263097429186324073070139212730,18781523847639554692492070299258169935566468,-121275574624841291294740106457630641456954640182250,0,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,-33046976527427969433650413380112451328120136191620,17325082089263041613534300922518663065279234311750,18008905658504805543025635497911412925086858851580,0,0,0,5402963,-112338294577728,-3081508888809517115520,82349590049879626345455133920,-255764122223332669757180975378478900,-5060051054516261255913825373053519712917936,0,1,-2371390,0,0,0,0,0,0,0,0,4644266,-144416129651496,-998103573133774880910,44373581533934249545163596960,-2806398470326794138888458631091170,-5366149670754158483569162942645191410161848,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,-15038070868923163890624777882201038403033277340040,0,-18008905658504805543025635497911412925086858851580,0,1,3973466,0,0,0,0,0,0,0,-5402963,112338294577728,3081508888809517115520,-82349590049879626345455133920,255764122223332669757180975378478900,5060051054516261255913825373053519712917936,0,0,0,-771100,139585311186948,-2427481790049864181440,-15921357666809978470936682680,716907274235989539559695275634122700,-2684559001397343842747054216222620003962972,0,739882,0,0,0,0,0,0,0,0,15433041046888888979327206349936258231366200572614,17325082089263041613534300922518663065279234311750],[0,0,-447426500232890640457842369370436667327162,0,843341842419376875985637562175586618819656,0,0,0,-4243028,0,0,0,0,0,0,0,-5402963,93615245481440,2054339259206344743680,-41174795024939813172727566960,85254707407777556585726991792826300,0,0,0,0,-771100,116321092655790,-1618321193366576120960,-7960678833404989235468341340,238969091411996513186565091878040900,0,1,8524144,0,0,0,0,0,-700958984478524379276471733774551146637668,894358278459026413928193823774198568360308,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[-18487311730143870741375508184190828452604723568352713050000041456,0,0,2233158801970430521035542921305517647815245309544838678371880412,21133520969921475892397488380709395343604653985475351794570724118,-771100,186113748249264,-4531299341426413138688,-44579801467067939718622711504,3345567279767951184611911286292572600,-25055884013041875865639172684744453370321072,-421065984329848588937493780701629075284931765521120,4513192727680294306047272477754397372222575422026923894800,-12225832,192554839535328,1863126669849713111032,-124246028295015898726458071488,13096526194858372648146140278425460,50084063593705479179978854131355119828177248,-485102298499365165178960425830522565827818560729000,-2030197176267001380200063703540941067670906163968646927232,1,5210176,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],[0,317969286500316286100423155589200855998127809164795016898,253774647033375172525007962942617633458863270496080865904,246179804459720502155485904130098815529694118588570469952,0,0,5402963,-131061343674016,-4314112444333323961728,144111782587289346104546484360,-596782951854442896100088942549784100,-17710178690806914395698388805687318995212776,126062339609533638801179448485379890475608011961060,1,-2371390,0,0,0,0,0,0,0,0,4644266,-168485484593412,-1397345002387284833274,77653767684384936704036294680,-6548263097429186324073070139212730,-18781523847639554692492070299258169935566468,121275574624841291294740106457630641456954640182250,0,0,0,0,0,0,0,0,0,0,0,0],[0,0,564149090960036788255909059719299671527821927753365486850,0,-246179804459720502155485904130098815529694118588570469952,1,3973466,0,0,0,0,0,0,0,-5402963,131061343674016,4314112444333323961728,-144111782587289346104546484360,596782951854442896100088942549784100,17710178690806914395698388805687318995212776,-126062339609533638801179448485379890475608011961060,0,0,-771100,162849529718106,-3398474506069809854016,-27862375916917462324139194690,1672783639883975592305955643146286300,-9395956504890703449614689756779170013870402,-105266496082462147234373445175407268821232941380280,739882,0,0,0,0,0,0,0,0,0,275991245103575021660367307310719176050140584341391469536,253774647033375172525007962942617633458863270496080865904],[0,0,-15038070868923163890624777882201038403033277340040,0,-18008905658504805543025635497911412925086858851580,0,0,-4243028,0,0,0,0,0,0,0,-5402963,112338294577728,3081508888809517115520,-82349590049879626345455133920,255764122223332669757180975378478900,5060051054516261255913825373053519712917936,0,0,0,-771100,139585311186948,-2427481790049864181440,-15921357666809978470936682680,716907274235989539559695275634122700,-2684559001397343842747054216222620003962972,1,8524144,0,0,0,0,0,0,15433041046888888979327206349936258231366200572614,-17325082089263041613534300922518663065279234311750,0,0],[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]):
x_theta_vars_filtered := [c, d, e, f, w_0, w_1, w_2, w_3, w_4, w_5, w_6, w_7]:
x_theta_vars := [a, b, x_0, z_0, y_0, z_8, z_7, z_6, z_5, z_4, z_3, z_2, z_1, y_8, y_7, y_6, y_5, y_4, y_3, y_2, y_1, x_9, x_8, x_7, x_6, x_5, x_4, x_3, x_2, x_1, w_7, w_6, w_5, w_4, w_3, w_2, w_1, w_0, f, e, d, c]:

# TODO: In below code, do not consider certain indices that correspond to loc. id. derivatives
#==============================================================================
GetBases := proc(Jac, r, ind_take, ind_left)
# we want to see all collections that are linearly independent among members 
# of Jac containing the elements of ind_take and maybe some elementes of ind_left
#==============================================================================
    local result, ind_left_new, ind_take_new, ind:

    if nops(ind_left) = 0 then
        if nops(ind_take) = r then
            return [ind_take]: # return the linearly independent columns 
        fi:
        if nops(ind_take) < r then
            return []:
        fi:
    fi:
    if nops(ind_take) + nops(ind_left) < r then
        return []:
    end if:
    if LinearAlgebra[Rank](Matrix(Jac[[op(ind_take), op(ind_left)]])) < r then
        return []:
    end if:
    ind := ind_left[1]:
    ind_left_new := ind_left[2..nops(ind_left)]:
    result := GetBases(Jac, r, ind_take, ind_left_new):
    ind_take_new := [op(ind_take), ind]:
    if LinearAlgebra[Rank](Matrix(Jac[ind_take_new])) = nops(ind_take_new) then
        result := [op(result), op(GetBases(Jac, r, ind_take_new, ind_left_new))]:
    end if:
    return result:
end proc:

#==============================================================================
GetTrBases := proc(Jac)
# returns the transcendence bases (as a list of indices of the parameters)
#==============================================================================
    local rankJacobian, JacCols, result, bases, b, complement, i:

    rankJacobian := LinearAlgebra[Rank](Jac):
    JacCols := [seq(Jac[.., i], i = 1..ArrayTools[Size](Jac)[2])]:
    bases := GetBases(JacCols, rankJacobian, [], [seq(i, i=1..nops(JacCols))]):
    result := []:
    for b in bases do
        complement := select(x -> not evalb(x in b), [seq(i, i = 1..nops(JacCols))]):
        result := [op(result), complement]:
    end do:
    return result:
end proc:

result := GetTrBases(JacX):

print(result);
