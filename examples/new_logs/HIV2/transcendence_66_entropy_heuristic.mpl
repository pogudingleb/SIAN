infolevel[Groebner]:=10;
Et_hat := [104648239370196012397-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 86145454481708475831-z_0, -c*q*w_0*y_0+h*z_0+z_1, 36303664780979656877-k_0, k_1, -w_1-259607321571041515488344681135811767763156894242874278787831934458592944702246283, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 364612464718126459599858196576707617587424647778553726378653565481775588494389607-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+644023844250977776388579573098734198726144775188870805191009186301496455723383982361233733736928238548808877276271642386453131645100689784165, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, -904516558965127558988926285815035165420015608706275917551732552696694893458344313855242973722228221187539609569679855092581338732665310343199-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3-1597669547429565906062707582549575531253188701974365332505953076399286304107507295619763234909912191447072792422916327713839331912292885979597836206550083653144683015256279730040302492585464176774762779, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 2243889840887936329948522542831213721005650137499094737300732194938671854220802724165750476362375804737552085494207040415658991662579007861280251631185927697262557426959106958841534323386201770681569953-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+3963437077011483006258924379880466469186319314261440604901207025676808229230505501640147084646098721727722092938108723700891022615848028657830886890601327226450002095824672548310525086187973271458825744525427284125979500012822077534047041346384828557832153456903, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, -5566555491035745203489009934042147079677229229776896323106678546128532589408840211624449043728950074987746253212736805719321649956245419843811713009321211907617551694860891252769893916835980878960533835689620349757918243413533375975670263955379299259860930510435-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5-9832342043886807185742627649946159039656071798771870179194688250907346975916797332990233881762921967606833161378203850183575109826955249066451908424807725123778139836481707894909521935704825715821700493256983629281576842303503540708224775484023059252061262872791698768866162475242778375981850669998238749234678433322698677, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 13809296459276463581337868264173048693515271629414064296872530870333547966877103359231097437773447544826989096004806161331297693644670768445604370727252329056431717145760967269144975862767127192596698161788277372776150359503507905706433240264157988121411669291120449800857779241270838802342105260824432310037392324919972237-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+24391695437455813811882803430971016369765645351034442056228061933008459832903538049902204789389210142129217214609963003030085909144974139027686435591128573326254105553443203720077482497126327556651513782121524786158423217788157413755392521266524956076110315811557258432498328057964696713868297936428121525629143603864806171121061211112643025995708622170264441278355709003317386636623, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7-60509978564416595990051725396934430260238019842472258482799350348519285254830788146641663492735579084774557347694230328653930523858911741188660498921319639669985958985303123510011622216305362527432142604691422460632803401518370482163271734573629003889238330023665477884685006172007125333081579324706662502178275743968991640551255107234454640656347549443688650031141749362206224251442872887030478133223130493732523300358772535793284849602406821, -34257570773753907046487939736860004076583292652803481092700542760574591561137865234926115542440807715939790636792498686336160110441706529429363727556131006484599311010913775994505305224928361038440120226948275684362048052990547944430346007711519375431188495118586190894006427575251541149853241786579534100579661037142222253781575782215004320190819917689825723974813998422249658146185-z_6, -k_1, -k_2, -k_3, -k_4, 691613346214019580335776435925576924298339383769600081446966314888514096444723444895783109345675827682328821137273952645022174614950357-y_3, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, y_3] 81
# 264, [1 = 9, -1 = 2, -15*c*w_2*x_4*y_0 = 1, c*q*w_3*y_2 = 1, -3*c*w_0*x_2*y_1 = 1, -164818727526407584224582241015472061561268250705468742805250854340282119696147721946275673955778194068305854889797040849026418884917373978479502188708810660198213701132365888216757381609938143300129045980020666701804658690708597812560949912633477767356876233539525765096819262942688165855214619303333905173103926951655537747511517048619957434971519104080160061107503551598056956786627727763745356294986395603386928265724434770341498707000415 = 1, -30*c*w_1*x_1*y_4 = 1, -4*c*w_1*x_0 = 1, c*q*w_2*y_1 = 1, w_7 = 1, -z_3 = 1, -6*c*w_5*x_1*y_0 = 1, -c*w_2*x_0*y_0 = 1, -5*c*w_4*x_0*y_1 = 1, a*y_1 = 1, z_aux = 1, -beta*v_0*x_5 = 1, beta*v_5*x_0 = 1, -10*c*w_0*x_2 = 1, -c*w_0*x_0*y_5 = 1, -4*q*w_1*c = 1, -c*w_6*x_0*y_0 = 1, -278915016615161777466941150036983400437558304405153015500624043150453025912641951665247663197050419117810744691810756103357604551176032186838 = 1, -y_0 = 1, w_3*q*c = 1, -beta*v_0*x_4 = 1, -20*c*w_1*x_3*y_1 = 1, -c*w_0*x_0 = 1, beta*v_3*x_0 = 1, w_1 = 1, -z_6 = 1, -2*c*w_1*x_1*y_0 = 1, -2*c*q*w_1*y_1 = 1, a*y_5 = 1, -5*c*w_1*x_4*y_0 = 1, -z_0 = 1, beta*v_0*x_5 = 1, c*q*w_3*y_1 = 1, beta*v_0*x_4 = 1, -4*c*w_0*x_1 = 1, beta*v_4*x_1 = 1, -z_4 = 1, x_5 = 1, beta*v_3*x_1 = 1, h*z_3 = 1, a*y_2 = 1, -10*c*w_3*x_2*y_0 = 1, w_3 = 1, -c*q*w_4*y_0 = 1, y_4 = 1, d*x_3 = 1, c*q*w_0 = 1, -c*q*w_0*y_1 = 1, d*x_0 = 1, -4*beta*v_1*x_3 = 1, beta*v_2*x_0 = 1, c*q*w_3*y_0 = 1, z_1 = 1, z_4 = 1, -2*c*w_1*x_0*y_1 = 1, beta*v_0*x_0 = 1, -60*c*w_2*x_3*y_1 = 1, -20*c*w_1*x_1 = 1, beta*v_0*x_3 = 1, -c*w_0*x_3*y_0 = 1, -c*q*w_0*y_2 = 1, -6*c*w_1*x_1*y_1 = 1, -3*beta*v_2*x_1 = 1, beta*v_1*x_1 = 1, x_6 = 1, -6*beta*v_2*x_2 = 1, -c*q*w_1*y_0 = 1, -w_6 = 1, -30*c*w_1*x_4*y_1 = 1, -3*c*w_2*x_1*y_0 = 1, -10*beta*v_2*x_3 = 1, h*z_5 = 1, -z_2 = 1, -5*c*w_1*x_0*y_4 = 1, u*v_4 = 1, beta*v_0*x_1 = 1, -c*w_4*x_0*y_0 = 1, -6*c*w_0*x_1*y_5 = 1, d*x_2 = 1, c*q*w_6*y_0 = 1, -10*beta*v_3*x_2 = 1, w_4 = 1, c*q*w_1*y_2 = 1, -15*c*w_4*x_2*y_0 = 1, c*q*w_2*y_4 = 1, u*v_2 = 1, y_5 = 1, -4*c*w_3*x_1*y_0 = 1, -90*c*w_2*x_2*y_2 = 1, z_5 = 1, -3*c*w_1*x_0*y_2 = 1, c*q*w_1*y_0 = 1, -60*c*w_1*x_3*y_2 = 1, beta*v_2*x_1 = 1, -c*w_3*x_0*y_0 = 1, -c*w_0*x_1*y_0 = 1, -c*w_0*x_2*y_0 = 1, c*q*w_0*y_0 = 1, -2*c*w_0*x_1*y_1 = 1, -5*beta*v_1*x_4 = 1, -c*q*w_5*y_0 = 1, -c*w_0*x_0*y_4 = 1, -w_3 = 1, -w_0 = 1, -6*c*w_5*x_0*y_1 = 1, c*q*w_4*y_0 = 1, -6*c*w_1*x_5*y_0 = 1, -y_1 = 1, -3*c*w_2*x_0*y_1 = 1, beta*v_1*x_3 = 1, -5*c*w_4*x_1*y_0 = 1, beta*v_3*x_2 = 1, -15*c*w_2*x_0*y_4 = 1, c*q*w_0*y_4 = 1, -3*c*w_1*x_2*y_0 = 1, -c*w_0*x_0*y_6 = 1, -4*beta*v_3*x_1 = 1, beta*v_4*x_0 = 1, -c*w_5*x_0*y_0 = 1, -w_5 = 1, -c*w_0*x_0*y_1 = 1, -beta*v_4*x_0 = 1, a*y_4 = 1, z_6 = 1, u*v_0 = 1, h*z_2 = 1, -4*c*w_3*x_0*y_1 = 1, -30*c*w_1*x_2*y_2 = 1, -c*w_0*x_0*y_0 = 1, -beta*v_3*x_0 = 1, -w_2 = 1, beta*v_1*x_2 = 1, a = 1, c*q*w_1*y_1 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_1*y_4 = 1, -60*c*w_1*x_2 = 1, -15*c*w_0*x_2*y_4 = 1, -beta*v_2*x_0 = 1, -c*q*w_0*y_4 = 1, -10*c*q*w_3*y_2 = 1, -6*c*w_0*x_2*y_2 = 1, -c*w_0*x_5*y_0 = 1, -beta*v_0*x_0 = 1, c*q*w_4*y_1 = 1, x_4 = 1, beta*v_0*x_2 = 1, c*q*w_5*y_1 = 1, -c*q*w_0 = 1, -10*c*w_2*x_3*y_0 = 1, -5*c*q*w_4*y_1 = 1, u*v_1 = 1, -c*w_0*x_0*y_2 = 1, c*q*w_1*y_5 = 1, y_6 = 1, -lm = 1, -15*c*w_4*x_0*y_2 = 1, -beta*v_5*x_0 = 1, -beta*x_0*v_1 = 1, h*z_1 = 1, -6*c*w_2*x_2*y_0 = 1, -60*c*w_2*x_1 = 1, -2*beta*v_1*x_1 = 1, -5*c*w_0*x_4*y_1 = 1, y_2 = 1, -200275591362100408122726486434978192082818297638257436524234088172441429289786695546364547781045103911678228724741854185416225592368885577976888230304127670336166398780833281003232908235933820396976353896937396568637104552052874068109358790706811413731355547934775775296348502446628354538461281317537642566573645375787349764862127053001630327821232612057741174355984636146724963318 = 1, -20*c*w_3*x_1*y_1 = 1, d*x_4 = 1, -4*c*w_0*x_3*y_1 = 1, -5*beta*v_4*x_1 = 1, q*w_2*c = 1, beta*v_2*x_2 = 1, -20*c*w_0*x_3 = 1, x_1 = 1, w_5 = 1, -6*c*w_2*x_0*y_2 = 1, -w_4 = 1, beta*v_2*x_3 = 1, -194503980712306719331673225396981017649589001127842111128364091688806229397220549852914261249712612573911644372217205957552089190787944608959325475955080188597633654975933405822389120551873648823906260156158992646533360385166232031578306153955993131239259904604743847353603901133386784210582875258498557781532618530608495 = 1, -beta*v_0*x_3 = 1, -12*c*w_2*x_1*y_1 = 1, -c*q*w_2*y_0 = 1, c*q*w_0*y_2 = 1, -15*c*w_0*x_4*y_2 = 1, c*q*w_0*y_5 = 1, -beta*v_0*x_2 = 1, -3*beta*v_1*x_2 = 1, z_3 = 1, -20*c*w_3*x_3*y_0 = 1, -10*c*w_2*x_0 = 1, b*w_2 = 1, -10*c*w_0*x_3*y_2 = 1, -c*w_1*x_0*y_0 = 1, b*w_0 = 1, -10*c*w_3*x_0*y_2 = 1, -60*c*w_3*x_2*y_1 = 1, b*w_3 = 1, q*w_1*c = 1, d*x_1 = 1, -6*c*w_0*x_5*y_1 = 1, z_2 = 1, -4*c*q*w_3*y_1 = 1, y_1 = 1, -30*c*w_2*x_1*y_2 = 1, -w_7 = 1, v_5 = 1, -10*q*w_2*c = 1, -4*c*w_1*x_3*y_0 = 1, x_2 = 1, -270877148049476526865509287601522425205673947385712555130824215342753506329726095 = 1, c*q*w_5*y_0 = 1, c*q*w_4*y_2 = 1, -5*c*q*w_1*y_4 = 1, -c*w_0*x_4*y_0 = 1, -c*q*w_0*y_0 = 1, b*w_6 = 1, c*q*w_2*y_0 = 1, c*q*w_0*y_6 = 1, -229535800213430833046071428647745561111868539496626746574080018173889401798672933180034287262201265441088228325350952475773216169046077712398889968292300765668791503946947052542571478451556870450681215 = 1, u*v_3 = 1, -z_1 = 1, -3*c*q*w_1*y_2 = 1, -z_5 = 1, w_6 = 1, v_4 = 1, a*y_0 = 1, -5*c*w_0*x_1*y_4 = 1, h*z_0 = 1, v_1 = 1, -236346926978904436719568454748778602012358613691752932997434020712583894484443859566116894969380093568709675867745002027093398462799396399447923208012074320225663661281767874529591683900145755749483687450867533603163647030756744060339181405835628716291944151078 = 1, h*z_4 = 1, -12*c*w_1*x_1*y_2 = 1, -w_1 = 1, -30*c*w_2*x_2*y_1 = 1, w_2 = 1, -c*w_0*x_6*y_0 = 1, b*w_1 = 1, beta*x_0*v_1 = 1, c*q*w_2*y_2 = 1, -c*q*w_3*y_0 = 1, b*w_4 = 1, beta*v_1*x_4 = 1, -60*c*w_3*x_1*y_2 = 1, -12*c*w_1*x_2*y_1 = 1, -6*c*w_1*x_0*y_5 = 1, d*x_5 = 1, v_3 = 1, -3*c*w_0*x_1*y_2 = 1, -y_4 = 1, v_2 = 1, b*w_5 = 1, -3*c*q*w_2*y_1 = 1, c*q*w_0*y_1 = 1, -6*c*q*w_2*y_2 = 1, -20*c*w_3*x_0 = 1, -beta*v_0*x_1 = 1, -y_2 = 1, x_3 = 1, -c*q*w_0*y_5 = 1]
# 273, -5.531957818
# 2
# [y_3 = 19, k = 5]
# [y_3 = [4, 4, 1, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4], k = [2, 2, 2, 2, 2]]