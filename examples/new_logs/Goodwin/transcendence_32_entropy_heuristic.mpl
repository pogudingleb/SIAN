infolevel[Groebner]:=10;
Et_hat := [4100027717355298791-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 2940692650853401964-gama_0, gama_1, -x1_1-141909603279363788161216643392101528938349695616846058751/4242963953987410304, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+1396905348258553439277461153152638605645218140508772671327076412904345801284898007199623350003060736741571257327/5119987319732837198871171193248909295179177755338582016, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-6875308319169655187354267960655529188470559408212551271782419316401767386455001705293591585515236845843658657926589816026739020797395859641044368728303347701447403807/3089145988335542775806755494694355696367404567331421210442938532368387856578366425837125632, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+67677975546829019121140064998155256616223210069215883292080816789017137779684675356776607749537358266713861815631326395651112473366129328103203685747121984414438244827565045978428978841527663438122204270797103241373374739/3727674441876362051724084961103474972623392916223077641350098533664876660413736670458276401832991752862946983954156643532668928, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-4193785209801254633424093368811535050309853643887487894534762049392237245855995806517247072953898988605072341295434133611542156529735021718623364500037302274110563499698280877790449528199384134134252300663030769777239756406002839923264851508466551783411375990253249932320539402348354149/2249093567783322438232097405470315677269886653461355485832816475656454957049065258027167702116199690128245810552027011581691029666937746768214198240533156154310656, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6+1307226685608318119451363630984000507499685301445474115355470413206482996091926191735746008368877675756147248453332638457060589175288883272987143784537330660222839679368511750954234099606800112174926502847579131832153944198856432868737665493744072122864145707158305661915843914113319754921326029847180023215614310723753108000722119294688174100468910279232109/1356991324086259898278146565714640652944789200774662896423938381113977499632584203729823629854944773087277076316148018592090752219215897438280036019218105166357562164329287021202451910431370391846912, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+12696520545248806036742432238614113201762871923935739301203153751018166039653178148040318499408335630941411677406402254994128114551753484235977237829941508553091987046155210664174646244470257256073053361772633023786830357997529979630292997115769510309063180047632416017077086676437700409291449663360615777730213878699522592350141405045271434198196648501617502610231413638938694485557842339513117496503759521115254314100677539177/102342643722258850860175571066128097332376396740246516853891694472764205075724510080372806520368203113525643011172392767710080759729562570838822597810827516610057735722923651223737123393475015175693086859701889917927423708433511088128, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8-265155941329561823559989813974969315556125288911543736097045065048389583256783233287133717065346290742954936224989840190522530729984957562294917759673097264219872311756134766448137400851641849926815356749971860967372353864974231410049159397741347644349793492901515594408460429002932364425300776145726442528014859378942630549599336916045122946741709487431913691172098771591026062624992555957653702611306907151337731332426117340904067071320237935039162505247322057415550375986714151055289159533548874281/493987734808322210547200839353826550267664918434588102048212627780839550348355488053422760823430852245462784766083616985244271896791187241473442255480111765922185265904820263024113858321769988922890506108644613431273905533812631334658853202987248338009871024768520552448, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9+89084115924955874895838219658978897514268592865459398814546675057859551580099607517744159034798189609745630044611590451966542553370422154525893043187986724068109206663263260316944433009381212942076308464938802254388690687193209395703403107574410464257115225053970944646319794381202498483472179516210046880298482166579898442636757099892344550762224070283874725553180402161747391299187662735896508286860296072008028457082732001638098730770229874231993459789463969546396015266209849329146161041716820600728161659029984235125589589716819123337288085466503168881712151945013141/298047657928519590286274943383965657861397152078494627373327659064734434959425699808245142884099734734232068868648591139833406496123494053276271536944324600772028495918787723573514983775575714088653907388159032843456007352240698276812717848097309725501506559383396368938456197194144877721574147973705629696, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, 507248149362102118971646063625064801330723050645767548627036415241081156236288405228508714838468653935626538191151782142770336055603327268987294865578853190061164679147138745611614959780863363640484920118982667645098679766408303864334075297323341283054497365/931918610469090512931021240275868743155848229055769410337524633416219165103434167614569100458247938215736745988539160883167232-x3_6, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x3_6] 173
# 275, [1 = 7, -1 = 2, x1_2*x4_0 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -3*sgm*x2_2*x4_1 = 1, b*x1_1*x4_2 = 1, -7*sgm*x2_1*x4_6 = 1, x3_1*x4_4 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, -5*sgm*x2_4*x4_1 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, -sgm*x2_0*x4_4 = 1, x1_5*x4_2 = 1, -5*sgm*x2_1*x4_4 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, -sgm*x2_4*x4_0 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, b*x1_1*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -sgm*x2_0*x4_6 = 1, x3_2*x4_3 = 1, -sgm*x2_0*x4_7 = 1, -sgm*x2_5*x4_0 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x4_2 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, delta*x4_0*sgm = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, -21*sgm*x2_2*x4_5 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -10*sgm*x2_3*x4_2 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, x1_7*x4_0 = 1, -x2_1 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -x2_6 = 1, -21*sgm*x2_5*x4_2 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, x1_3*x4_5 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, delta = 1, x3_4 = 1, -sgm*x2_3*x4_0 = 1, b*x1_1*x4_5 = 1, -3*sgm*x2_1*x4_2 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, x1_8*c = 1, x1_2*x4_3 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, -20*sgm*x2_3*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, -x2_0*x4_0*sgm = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, -x2_2 = 1, -sgm*x2_0*x4_3 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, -4*sgm*x2_3*x4_1 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, -35*sgm*x2_3*x4_4 = 1, x1_3*c = 1, -sgm*x2_1*x4_0 = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, -sgm*x2_0*x4_2 = 1, x1_1*x4_8 = 1, -x2_5 = 1, x2_4 = 1, -35*sgm*x2_4*x4_3 = 1, x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, -6*sgm*x2_1*x4_5 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -x1_1 = 1, -x2_0 = 1, -x1_6 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, delta*x4_1*sgm = 1, -x2_3 = 1, x3_1*x4_6 = 1, -sgm*x2_6*x4_0 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, -sgm*x2_0*x4_1 = 1, -sgm*x2_7*x4_0 = 1, x3_0*x4_4 = 1, x3_1*x4_5 = 1, x3_3*x4_5 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -x1_4 = 1, x2_1 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -6*sgm*x2_5*x4_1 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, -2*sgm*x2_1*x4_1 = 1, b*x1_7*x4_1 = 1, b*c*x1_4 = 1, -4*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -6*sgm*x2_2*x4_2 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -x2_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, x3_4*x4_1 = 1, -10*sgm*x2_2*x4_3 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, -15*sgm*x2_4*x4_2 = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, -sgm*x2_0*x4_5 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x3_7*x4_1 = 1, x3_1*x4_3 = 1, x3_3*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, -7*sgm*x2_6*x4_1 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, -sgm*x2_2*x4_0 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.588688399
# 2
# [gama = 43, x3_6 = 6]
# [gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2], x3_6 = [4, 2, 1, 4, 2, 2]]