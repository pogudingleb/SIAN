infolevel[Groebner]:=10;
Et_hat := [6919518019486802806-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-120246482962306050998535079380136468129404759385584165999/5976389535620664620, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+2908666204922832242397915398847849495835811567399697822501096444678371833609257058913096557249811221863064754829/49716753035278232092388475610811528335845012165073019500, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-70358308062218119125080214758522730040496871009571805481102455861541518149187150462455989969178260034909052483321978926657271141823779213662167613373140873068742441383/413586751271585017233560039694425463241164311326554950450629520761044650762065804990632387500, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0)*gama)*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+3428908911197872886804481149705741356335875550266922712932024207151845200086384574776559518834277902102246869035490528487535250686731362992726984475664828414058040902837764144354375099646072917707934903120914116823628767273/6881141280727109483138189106387256748065993899488804965668601161126545245874027587155214001149701630169612202995845466049924375000, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+94255474333074764449316449409944365334017977982587731720297913282991213144211442239901125489566521080993065492506464118795496857844355826284131844632522122746794868524439259440898508318909473947699323922549616989454383100088686672398460411446996079924995888493029597771590219182747442993020883/57243256922214491735595450275905343493968564413801569752337714453930737071267262551695475565111943427835240116038125574048607066565824354310872999379514744836734375000, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-3794035785204290430228114522657096173847413281388424360040218434738213794955586140242421017923773408539237225570274483968429000588461357259822607839110993187684033913291001963587617416571767614655593761282441605428578276798644182718735659504807617963856045798834821098435500599958638131052814445322906080517108168181527734300939587674894506204451769075639266703261526/11904967248094376947014083231376116636616908833541230970515922270566239632764358765388554601293763335088185201266315960831048122046995726057815200726092250064119138343035389882591891502623996915224609375, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+6868287026677964042542324357307392968171135486118421605101053413458747118694554745634854459719537040009953453187031093347606569492275606851479880669309664148161649029369572753240492535789547102470069013942105067510048865015564218135053605459020843869253858785283193911191402345105312146009241084771239975191431490807113984571871219395089159199505044151323950362625142324547052027044239879370608551070822774172560772071895473636189423050717371/1584572258655803461299189641217892315103278519816290759568937139419487518481332765102386993997844084463772090516492620840490968183975005214822067570306519647919633257464448336568469919260547326383974516446866002094303952700921820746093750000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-222575373548814349419758383602119034871574732549659567940918444089742780322611792534438846199395412041681015992577477774911608271625199087911555783063498636179068082526822044946713842938088530377682716247562957103116591018951584276920710684471262706272964481921854699557590257872872677414463793888417233687098154352909959531514014811930432027037281695733101404766754316828324658741051917816972754977579642029710970878594785274966826764633096068034978347907568148485496395249285322881793857183715076482405520654045747/13181836153851390692966952609820599558463018098179861705250666119057654788645680800223338721218455220511157124473953907223404295977105999906230697176396802728388540811098263693328224327828812534102361540649633863399076361626541417864346049246222383223713112630037378808593750000, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-73632116493735916352244443091234420248406059273891445513734191386836784311537449904553018708919497035189189765243845795682810439360493241753317198401417187956779207670359813701850722489765898399293359028334383854897420068784267237106407161980652536262032771043668485991620684200319281530019665212090683602179813969664440545565513727269664172900786992726209151383165350919793188429439297038342440914059504600713936207090765945634201438079507033031468644725106061738416408295621084865295356079435966042445097917245194748407436523395144203239252662180620422087107486697340806144410841993089449/109657860938689751995661859788698693879079704173981116374454605390898764483301132478280123130893351171847571391132127947284958999183368043295475035252271212312224549399133121282143192340196806738773535981950526418642001238278585634841276603774372407769947600897151249753059048696003746795968112245254567871093750000, -354565899167234961943045790578069740349663258505585446010218083029210451992173648450369886761340910583974230919951939745329645817677413749676600244589760083372967850226314812210517899403/206793375635792508616780019847212731620582155663277475225314760380522325381032902495316193750-x2_4, 5567817035794706995-x3_0, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_4, x3_0] 77
# 275, [1 = 7, -1 = 2, x3_1*x4_5 = 1, delta*x3_2 = 1, x2_2 = 1, x3_1 = 1, x2_5 = 1, x1_3*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_5 = 1, -gama = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, -gama*sgm*x2_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, -gama*x2_2 = 1, -21*gama*sgm*x2_5*x4_2 = 1, x4_7*delta*sgm = 1, x3_4 = 1, delta*sgm*x3_4*x4_1 = 1, x3_1*x4_4 = 1, b*c*x1_7 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_1*x4_3 = 1, x3_5*x4_2 = 1, x4_6 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, b*x1_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_7*x4_0 = 1, x4_4*delta*sgm = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -gama*x2_6 = 1, -35*gama*x4_3*sgm = 1, x1_3*x4_5 = 1, -gama*sgm*x2_0*x4_3 = 1, -alpha*x1_2 = 1, x1_6*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, -x1_4 = 1, -gama*x2_0 = 1, x1_1*c = 1, x1_1*x4_6 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, x4_3*delta*sgm = 1, delta*x4_0*sgm = 1, -gama*x2_5 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_0*x4_6 = 1, x1_2*c = 1, x3_7 = 1, c*z_aux = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, -gama*sgm*x2_7*x4_0 = 1, x4_7 = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, delta*sgm*x3_4*x4_3 = 1, x4_1 = 1, x3_1*x4_7 = 1, x1_1*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, b*x1_3*x4_0 = 1, beta*x2_5 = 1, delta*sgm*x3_3*x4_2 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, -4*gama*sgm*x2_3*x4_1 = 1, delta*x3_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, x3_3*x4_2 = 1, x1_3*x4_4 = 1, delta*sgm*x3_1*x4_3 = 1, x1_1*x4_4 = 1, -x1_0 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, x2_3 = 1, x1_9*x4_0 = 1, -15*gama*sgm*x2_2*x4_4 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, x2_1 = 1, delta*sgm*x3_4*x4_0 = 1, x4_2*delta*sgm = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, -20*gama*sgm*x2_3*x4_3 = 1, x3_2*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, x4_0*z_aux = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_2 = 1, delta = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*x3_5 = 1, x1_5*x4_4 = 1, -gama*sgm*x2_0*x4_7 = 1, x1_5*x4_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_0*x4_2 = 1, -alpha*x1_5 = 1, delta*sgm*x3_5*x4_1 = 1, x3_3*x4_5 = 1, b*c*x1_5 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -gama*sgm*x4_0 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x3_1*x4_6 = 1, b*x1_0*x4_7 = 1, -5*gama*x4_1*sgm = 1, b*x1_1*x4_2 = 1, x1_3*c = 1, x1_4*x4_2 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x1_5*x4_0 = 1, x4_5*delta*sgm = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x4_6*delta*sgm = 1, x3_3 = 1, beta*x2_0 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, beta = 1, b*x1_5*x4_1 = 1, x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_1*x4_7 = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, x4_4 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, x1_3*x4_0 = 1, b*x1_1*x4_6 = 1, beta*x2_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_1 = 1, b*x1_1*x4_4 = 1, x3_1*x4_1 = 1, x3_3*x4_1 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_5*c = 1, -3*gama*sgm*x2_1*x4_2 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_6*x4_2 = 1, -10*gama*sgm*x2_3*x4_2 = 1, b*x1_1*x4_5 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*sgm*x3_1*x4_6 = 1, b*c*x1_1 = 1, x3_5 = 1, b*c*x1_6 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, x4_3 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x1_8*x4_1 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, b*x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, -alpha*x1_3 = 1, x3_7*x4_1 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, -15*x4_2*gama*sgm = 1, delta*sgm*x3_2*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_4*x4_0 = 1, -gama*sgm*x2_2*x4_0 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_4*x4_4 = 1, x1_4*x4_0 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, delta*sgm*x3_5*x4_0 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x4_5 = 1, -gama*x2_1 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, x4_8 = 1, x3_5*x4_3 = 1, x1_7*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, delta*x4_1*sgm = 1, -gama*x2_3 = 1, x3_6*x4_2 = 1, x3_2 = 1]
# 282, -5.588688399
# 2
# [x3_0 = 19, x2_4 = 7]
# [x3_0 = [4, 2, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 3, 3], x2_4 = [4, 1, 4, 2, 2, 4, 4]]