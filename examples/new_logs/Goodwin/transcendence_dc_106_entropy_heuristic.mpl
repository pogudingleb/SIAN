infolevel[Groebner]:=10;
Et_hat := [3889358523717833500-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-173104665288220680512420881085948879666000133176405501999/5034329578474006803, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+3118814742670339947317133198895204625549332904126327950512331958892390151710952534204539970417632701485431231514/10259667812237971510045526793933336512598740814046823757, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-168574406407285131343677173627875302851917960952097252369246408713885200905974956942746636762174177526876898045703051827590597892916163375037966986721174318803157030323/62725800115004548060353927939844602144674397213765864451812337025636935926962135214418790649, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0)*gama)*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+337462136493688825283446163165621716321428729768920682649217427276974189714203385937419457673388704415699730323643220563106480434423784852094010527518408189228739167838855102928378188211491023837971335325624908983783276170/14203499358318729786505730641102269366528332575162970454668634516350811218483761695855775993130688375311126898229568710694851759, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+69122963668759521459963503291981212970006605530405579431458464382840347070072449934046046503985527997301704495540421641165350842075597058891152307335544246453152675369303764591075861297394677868141509279615030236203929175455778404355203132725474870361294816798439942367697256088893888354897/521026145088700914170327552769117222652346046474131912056972468635126922920082004153299217328814155333140466377105556680817294193755818431855841408973145759052254578, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+31115013882036952588090419464179032831216301113481509176998153535608920604242874301910757379338760983288744505173825670791256718524765640442083253964343568719341547543479534422143920252402808163461263917979655396406746380724782656495589638368539843102084643400821936212201529244542246554204063563460638487617108032760305011525754007367385265466711287670481922295/18307252808033783299721655819612057526388308899457139137687530080936364997385848944083070801397582978978426196357141746286521320070097879562082192642742174572232730170256651413785251615929917297559379, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-11780791916324934387846308786958868724637715802024722582604325853525505315687440587475150989996595253075354693825263151378200182265816577953277009804804095389536337039515570859783976189710235508545881828803393856387019708839819848684081152890452279851396829795326642142145397345932707482326153557339912471405623837777975969546554369637482919517065086356662545258679854847035513203966810903613475425239557989509724842070827416605246097559/6491784322438347562070369615715505282472009128021550648275280729929927205830391579429261167987780064510732478476716571561325431607491839986464759900570490020091983399447174826716043470618358199049934057758608183091873333404665462532405974, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-150764205810635178156058964514397305120704776191434540146171538227504556927755238510881379767435842586034136050898984993036501928665960620969389685283389444954860766914152656966246754794769134280805552514282615136414885084333645639347254838360659638045951828319201044594884524209909327278427197802959645857301544868514588482599326662907008654111717082887436557822943057510212889528972415792705932620306310376080180578059361459543065909133431228595210876533814766118469402644306812874343890629010073415585616346/6614937462745585946381661408386739240250874006432756398411557212703766213075270050034858614369587918812835472020772020103302011204165133561638156943165272759427765370954203380479256970640747466146411479037945586608703730760016702724300722153856973272937223708649016577547253, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+1718561383969286662146543599519117222934470038228848315487512076320981794751190483657795242653677767792041847187130542501174821135715330413103788118285451846943512354161622991785249722231048227401393440221237696991384503809148004729901772787880597429112218119568761745510810500436064819805245317932841983172354895376387778106783408768991476536409483829693458273068911104714228825976519162447541178086602011920217349916722269114867079397318718484243384840474697864686323914661718118521246282065623348601793525894081400521815124975368921441439179726241432080857591533365933915894618533/40442561363097941100124766047667420025013412756176351318073960445474919412015558663572543169133419201778816773073942421687863992837190449870736670406819703203784473875579843869106308275629938002608799555023356903564414292004534847030567055301906288566938668400797701483044044344188083302408427037818248494801921, 5310174902707922722585217352844028679234748845106217599602889395793675855001471689831048400039948568333901770034437481381439463710155903291959763244784546521420355116371576160851362822471978667525381611273492332364238569957289104206178872/165799603404498791282168061958392249414727617608906269898855655832888846908370759095592715873120875976394477411240101681262861-x2_5, -1191021010640278303780558405095929993244901287424428888452729200313124549556692408330419928067602488310190246661730412449038129980186238784735362081674482869218918918509510858478583273294740552689954994/81356420382625872970627662697593517697372759032121743776669697828322874094633119603656019-x3_5, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_5, x3_5] 38
# 274, [1 = 8, -1 = 2, x3_1*x4_5 = 1, delta*x3_2 = 1, x2_2 = 1, x3_1 = 1, x1_3*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_5 = 1, -gama = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, delta*x3_0 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, delta*sgm*x3_0*x4_0 = 1, -gama*sgm*x2_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x3_0*x4_3 = 1, -gama*x2_2 = 1, z_aux*x3_0*x4_0 = 1, x3_4 = 1, delta*sgm*x3_4*x4_1 = 1, x3_1*x4_4 = 1, b*c*x1_7 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_1*x4_3 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, x3_0*x4_5 = 1, b*x1_6*x4_2 = 1, b*x1_7*x4_0 = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -gama*x2_6 = 1, x1_3*x4_5 = 1, -gama*sgm*x2_0*x4_3 = 1, -alpha*x1_2 = 1, x1_6*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, -x1_4 = 1, -gama*x2_0 = 1, x1_1*c = 1, x1_1*x4_6 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, x4_2*delta*sgm = 1, delta*x4_0*sgm = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_0*x4_6 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, x3_7 = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, -gama*sgm*x2_7*x4_0 = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, delta*sgm*x3_4*x4_3 = 1, x4_1 = 1, x3_1*x4_7 = 1, delta*sgm*x3_0*x4_7 = 1, x1_1*x4_7 = 1, x3_0*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, b*x1_3*x4_0 = 1, delta*sgm*x3_3*x4_2 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, -4*gama*sgm*x2_3*x4_1 = 1, delta*x3_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_3*x4_2 = 1, x1_3*x4_4 = 1, delta*sgm*x3_1*x4_3 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_4*x4_0 = 1, -x1_0 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, x2_3 = 1, x1_9*x4_0 = 1, -15*gama*sgm*x2_2*x4_4 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, x2_1 = 1, delta*sgm*x3_4*x4_0 = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, -20*gama*sgm*x2_3*x4_3 = 1, x3_2*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_2 = 1, x2_4 = 1, -6*gama*sgm*x2_1*x4_5 = 1, x1_5*x4_4 = 1, -gama*sgm*x2_0*x4_7 = 1, x1_5*x4_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_0*x4_2 = 1, -alpha*x1_5 = 1, x3_3*x4_5 = 1, b*c*x1_5 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -gama*sgm*x4_0 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x3_1*x4_6 = 1, b*x1_0*x4_7 = 1, b*x1_1*x4_2 = 1, x3_0*x4_6 = 1, x1_3*c = 1, x1_4*x4_2 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x1_5*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x3_3 = 1, beta*x2_0 = 1, x3_0*x4_8 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, x3_0*x4_1 = 1, z_aux*x3_0*c = 1, beta = 1, b*x1_5*x4_1 = 1, x4_2 = 1, x4_1*delta*sgm = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_1*x4_7 = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, x1_3*x4_0 = 1, b*x1_1*x4_6 = 1, beta*x2_1 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_1*x4_1 = 1, b*x1_1*x4_4 = 1, x3_1*x4_1 = 1, x3_3*x4_1 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_0*x4_6 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_5*c = 1, -3*gama*sgm*x2_1*x4_2 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_6*x4_2 = 1, -10*gama*sgm*x2_3*x4_2 = 1, b*x1_1*x4_5 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*sgm*x3_1*x4_6 = 1, b*c*x1_1 = 1, b*c*x1_6 = 1, x4_3 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x1_8*x4_1 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, b*x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, -alpha*x1_3 = 1, x3_7*x4_1 = 1, -gama*x2_4 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, -6*gama*x4_1*sgm = 1, delta*sgm*x3_2*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_4*x4_0 = 1, -gama*sgm*x2_2*x4_0 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_5 = 1, x1_4*x4_0 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, x3_0*x4_4 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -35*gama*sgm*x2_3*x4_4 = 1, -gama*x2_1 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, -21*x4_2*gama*sgm = 1, x1_7*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, x3_0*x4_2 = 1, -gama*x2_3 = 1, x3_6*x4_2 = 1, x3_2 = 1]
# 282, -5.577999884
# 2
# [x3_5 = 8, x2_5 = 6]
# [x3_5 = [4, 2, 1, 4, 2, 2, 4, 2], x2_5 = [4, 1, 4, 2, 2, 4]]