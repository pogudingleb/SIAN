infolevel[Groebner]:=10;
Et_hat := [44076786280489047518-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 257249975016498225427-z_0, -c*q*w_0*y_0+h*z_0+z_1, -w_1+198678477824876538653303655640077047449440674896274306016245210949161632162206, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 737649914585580397810080585162320029986537262773995029227795653755878238414197-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+895553893144952774969199634230520204346142911056483951239871167123691596049867644387760701797085838813893928827585864345494752288697982, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k*y_0+u*v_0+v_1, 3324996547272942893954378182589228402971501223904911029226227175735106509250718286053418966745011287382530244097986913580896027004018755-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+4036757198401798762656863421218533722208783371877125513390425630417593017688565068351793117789740929225657577057956410302050960744255556202956904942544293862788751918879740710212772324480718046, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k*y_1+u*v_1+v_2, 14987600243386676215943256696369418800198052375593920899710198526861184426706215631461515133233033352493167275214427665706279079822490321715787804124503326173427734291841383401635526284914747781-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+18195899547288567236785441808792091173177740982120487917542318349070968134496633062367244595917601785736949206878035797341702838098486871563584977925545993783259080019023155230158181857350028526802862209548704905867010453338472657030146288723110689470, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k*y_2+u*v_2+v_3, 67557411823419053797535683126673467790743713415061168236957148176668798952022238764224683163042189064741196280825896101058137188294429215579236889300388574104903816091834899538390951582257496156595731628959660015368324132682927798619764310612377884083-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5+82018992984294315533377102730591373246898070371237097593352799601604353946949720989850488528950723890564302919013211288200153591667566620018696815384131424779889123295166282547224402853469562629956585870206178953455508748501167713847461358144364239581945550816826464435291641071358221193275624915942625056286, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k*y_3+u*v_3+v_4, 304518656633701373792957335255954276240670561172669776441598216433662362672967238099484799958955459901996704435402739615547422720834757685043211905967657122646019693341070549480687690143938696821905116430589932633908675651282492849589515161790818120947033469581940195644067485210058530165008071044075446094485-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+369705009234355716152054917963885073435059831249293279591126273956974068515617354526453572594412230674058156333491677873494465759639141124037846783300950266702062764123960526576666721379145189227255524827228310653178350265845514601885282884514674532499731989221730838139724285972917773613021233910543546159863441726742508223290610217854359699635923171865405033251070, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k*y_4+u*v_4+v_5, -w_7+1666465154956830497287342274295560128229738330362872157109903333089912213354717603369235258374425508022483173159729538306949859870892677083012036370354578623447005277798662367383015131916351999766323309291090398065960832436988400689449688067117108683045589317927624540605062865579584510369109254431730114037057832519816823459570851547601552701848994563685130253287504801032396213404497040389607720694874987704210347682095454, 1372634174920366658904799719174720607733456471619642637946916633330876839651659845134335224817853766668760767644872871313762131906309656239023856398631536218711645074527099728346377442943644341522748884512764489258897252883708749625115174748413759130538491860684679399689154829170048867566336474729678634243021795169135775223602701186893482928424184206893944402118179-z_6, -29121454681162707102135124390390696343376-v_1, 1997738662804165644-y_0, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, z_3, y_3, x_3, w_3, v_3, z_2, y_2, x_2, w_2, v_2, z_1, y_1, x_1, w_1, v_1, z_0, y_0, x_0, w_0, v_0, z_aux, w_aux, a, b, beta, c, d, h, k, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [v_1, y_0] 201
# 265, [1 = 9, -c*w_0*x_5 = 1, c*q*w_3*y_2 = 1, -3*c*w_0*x_2*y_1 = 1, -164818727526407584224582241015472061561268250705468742805250854340282119696147721946275673955778194068305854889797040849026418884917373978479502188708810660198213701132365888216757381609938143300129045980020666701804658690708597812560949912633477767356876233539525765096819262942688165855214619303333905173103926951655537747511517048619957434971519104080160061107503551598056956786627727763745356294986395603386928265724434770341498707000415 = 1, -30*c*w_1*x_1*y_4 = 1, c*q*w_2*y_1 = 1, w_7 = 1, -z_3 = 1, -5*c*w_4*x_0*y_1 = 1, a*y_1 = 1, z_aux = 1, -beta*v_0*x_5 = 1, beta*v_5*x_0 = 1, -c*w_0*x_0*y_5 = 1, -c*w_4*x_0 = 1, c*q*w_6 = 1, c*q*w_3*y_3 = 1, -278915016615161777466941150036983400437558304405153015500624043150453025912641951665247663197050419117810744691810756103357604551176032186838 = 1, -beta*v_0*x_4 = 1, -20*c*w_1*x_3*y_1 = 1, -k*y_2 = 1, beta*v_3*x_0 = 1, q*w_1*c = 1, w_1 = 1, -z_6 = 1, -6*c*w_2*x_2 = 1, -2*c*q*w_1*y_1 = 1, a*y_5 = 1, -z_0 = 1, beta*v_0*x_5 = 1, c*q*w_3*y_1 = 1, beta*v_0*x_4 = 1, beta*v_4*x_1 = 1, u = 1, -z_4 = 1, x_5 = 1, beta*v_3*x_1 = 1, beta*x_2 = 1, h*z_3 = 1, a*y_2 = 1, w_3 = 1, -c*w_0*x_6 = 1, y_4 = 1, d*x_3 = 1, -c*q*w_0*y_1 = 1, d*x_0 = 1, beta*v_2*x_0 = 1, -60*c*w_1*x_2*y_3 = 1, -4*c*w_3*x_1 = 1, z_1 = 1, z_4 = 1, -2*c*w_1*x_0*y_1 = 1, -c*q*w_0*y_3 = 1, beta*v_0*x_0 = 1, -60*c*w_2*x_3*y_1 = 1, beta*v_0*x_3 = 1, -5*c*w_1*x_4 = 1, -c*q*w_0*y_2 = 1, -60*c*w_2*x_1*y_3 = 1, -6*c*w_1*x_1*y_1 = 1, -3*beta*v_2*x_1 = 1, x_6 = 1, -6*beta*v_2*x_2 = 1, -c*w_0*x_3 = 1, -w_6 = 1, -30*c*w_1*x_4*y_1 = 1, -4*c*w_1*x_3 = 1, -10*beta*v_2*x_3 = 1, h*z_5 = 1, -z_2 = 1, -c*q*w_0 = 1, -5*c*w_1*x_0*y_4 = 1, u*v_4 = 1, -6*c*w_5*x_1 = 1, beta*v_0*x_1 = 1, -15*c*w_2*x_4 = 1, -q*w_5*c = 1, -k = 1, -c*w_6*x_0 = 1, -6*c*w_0*x_1*y_5 = 1, -k*y_1 = 1, d*x_2 = 1, -10*beta*v_3*x_2 = 1, w_4 = 1, -3*beta*x_2 = 1, c*q*w_1*y_2 = 1, c*q*w_2*y_4 = 1, u*v_2 = 1, -6*c*w_1*x_5 = 1, y_5 = 1, -3*c*w_1*x_2 = 1, w_3*q*c = 1, q*w_2*c = 1, -90*c*w_2*x_2*y_2 = 1, z_5 = 1, -3*c*w_1*x_0*y_2 = 1, -60*c*w_1*x_3*y_2 = 1, beta*v_2*x_1 = 1, -20*c*w_3*x_3 = 1, -2*c*w_0*x_1*y_1 = 1, -15*c*w_4*x_2 = 1, -c*w_0*x_0*y_4 = 1, -w_3 = 1, -w_0 = 1, -20*c*w_1*x_1*y_3 = 1, q*w_4*c = 1, -6*c*w_5*x_0*y_1 = 1, -3*c*w_2*x_0*y_1 = 1, beta*x_0 = 1, beta*v_3*x_2 = 1, -15*c*w_2*x_0*y_4 = 1, c*q*w_0*y_4 = 1, -c*w_0*x_1 = 1, -c*w_0*x_0*y_6 = 1, -4*beta*v_3*x_1 = 1, beta*v_4*x_0 = 1, -w_5 = 1, -c*w_0*x_0*y_1 = 1, -beta*v_4*x_0 = 1, -3*c*w_2*x_1 = 1, a*y_4 = 1, z_6 = 1, u*v_0 = 1, -c*w_0*x_4 = 1, h*z_2 = 1, -beta*x_0 = 1, -4*c*w_3*x_0*y_1 = 1, x_4*beta = 1, -30*c*w_1*x_2*y_2 = 1, -beta*v_3*x_0 = 1, -w_2 = 1, c*q*w_1*y_1 = 1, a = 1, -k*y_3 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_1*y_4 = 1, -q*w_1*c = 1, -15*c*w_0*x_2*y_4 = 1, -beta*v_2*x_0 = 1, c*q*w_0*y_3 = 1, -w_3*q*c = 1, -c*q*w_0*y_4 = 1, -10*c*q*w_3*y_2 = 1, -6*c*w_0*x_2*y_2 = 1, -c*w_0*x_2 = 1, c*q*w_0 = 1, -beta*v_0*x_0 = 1, c*q*w_4*y_1 = 1, x_4 = 1, beta*v_0*x_2 = 1, c*q*w_5*y_1 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, -c*w_3*x_0 = 1, -c*w_0*x_0*y_2 = 1, c*q*w_1*y_5 = 1, -q*w_2*c = 1, y_6 = 1, -lm = 1, -15*c*w_4*x_0*y_2 = 1, -beta*v_5*x_0 = 1, h*z_1 = 1, -4*c*q*w_1*y_3 = 1, -20*c*w_0*x_3*y_3 = 1, -5*c*w_0*x_4*y_1 = 1, -q*w_4*c = 1, y_2 = 1, -200275591362100408122726486434978192082818297638257436524234088172441429289786695546364547781045103911678228724741854185416225592368885577976888230304127670336166398780833281003232908235933820396976353896937396568637104552052874068109358790706811413731355547934775775296348502446628354538461281317537642566573645375787349764862127053001630327821232612057741174355984636146724963318 = 1, -20*c*w_3*x_1*y_1 = 1, -10*c*q*w_2*y_3 = 1, d*x_4 = 1, -4*c*w_0*x_3*y_1 = 1, -5*beta*v_4*x_1 = 1, a*y_3 = 1, -2*c*w_1*x_1 = 1, beta*v_2*x_2 = 1, c*q*w_1*y_3 = 1, x_1 = 1, -c*w_5*x_0 = 1, w_5 = 1, -6*c*w_2*x_0*y_2 = 1, -w_4 = 1, beta*v_2*x_3 = 1, -194503980712306719331673225396981017649589001127842111128364091688806229397220549852914261249712612573911644372217205957552089190787944608959325475955080188597633654975933405822389120551873648823906260156158992646533360385166232031578306153955993131239259904604743847353603901133386784210582875258498557781532618530608495 = 1, -beta*v_0*x_3 = 1, -12*c*w_2*x_1*y_1 = 1, -4*beta*x_3 = 1, c*q*w_0*y_2 = 1, -15*c*w_0*x_4*y_2 = 1, c*q*w_0*y_5 = 1, -beta*v_0*x_2 = 1, z_3 = 1, -c*w_1*x_0 = 1, b*w_2 = 1, -10*c*w_0*x_3*y_2 = 1, b*w_0 = 1, -c*w_2*x_0 = 1, -10*c*w_3*x_0*y_2 = 1, -2*beta*x_1 = 1, -60*c*w_3*x_2*y_1 = 1, b*w_3 = 1, d*x_1 = 1, -6*c*w_0*x_5*y_1 = 1, z_2 = 1, -4*c*q*w_3*y_1 = 1, y_1 = 1, -30*c*w_2*x_1*y_2 = 1, -w_7 = 1, v_5 = 1, x_2 = 1, -270877148049476526865509287601522425205673947385712555130824215342753506329726095 = 1, c*q*w_4*y_2 = 1, c*q*w_2*y_3 = 1, -c*w_0*x_0 = 1, -5*c*q*w_1*y_4 = 1, b*w_6 = 1, -10*c*w_2*x_0*y_3 = 1, c*q*w_0*y_6 = 1, -229535800213430833046071428647745561111868539496626746574080018173889401798672933180034287262201265441088228325350952475773216169046077712398889968292300765668791503946947052542571478451556870450681215 = 1, u*v_3 = 1, -z_1 = 1, beta*x_1 = 1, -3*c*q*w_1*y_2 = 1, -z_5 = 1, -5*c*w_4*x_1 = 1, w_6 = 1, v_4 = 1, -10*c*w_0*x_2*y_3 = 1, beta*x_3 = 1, -5*c*w_0*x_1*y_4 = 1, h*z_0 = 1, q*w_5*c = 1, -236346926978904436719568454748778602012358613691752932997434020712583894484443859566116894969380093568709675867745002027093398462799396399447923208012074320225663661281767874529591683900145755749483687450867533603163647030756744060339181405835628716291944151078 = 1, h*z_4 = 1, -12*c*w_1*x_1*y_2 = 1, -w_1 = 1, -30*c*w_2*x_2*y_1 = 1, w_2 = 1, b*w_1 = 1, -10*c*w_2*x_3 = 1, c*q*w_2*y_2 = 1, b*w_4 = 1, -60*c*w_3*x_1*y_2 = 1, -c*w_0*x_0*y_3 = 1, -10*c*w_3*x_2 = 1, -12*c*w_1*x_2*y_1 = 1, -6*c*w_1*x_0*y_5 = 1, d*x_5 = 1, -20*c*w_3*x_0*y_3 = 1, v_3 = 1, -3*c*w_0*x_1*y_2 = 1, v_2 = 1, b*w_5 = 1, -5*x_4*beta = 1, -3*c*q*w_2*y_1 = 1, y_3 = 1, c*q*w_0*y_1 = 1, -6*c*q*w_2*y_2 = 1, -k*y_4 = 1, -beta*v_0*x_1 = 1, x_3 = 1, -1 = 1, -c*q*w_0*y_5 = 1, -4*c*w_1*x_0*y_3 = 1]
# 273, -5.537035820
# 1
# [v_1 = 12, y_0 = 43]
# [v_1 = [3, 3, 1, 3, 3, 2, 3, 3, 3, 3, 3, 3], y_0 = [4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]]