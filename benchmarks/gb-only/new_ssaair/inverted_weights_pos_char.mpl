kernelopts(printbytes=false, assertlevel=1):
interface(echo=0, prettyprint=0):
infolevel[Groebner]:=10;
et_hat:=[24755770674423839499359380-Sd_0, Ad_0^2*Sd_0*b_a^2*eps_a^2*eps_s^2+An_0^2*Sd_0*b_a^2*eps_s^2+In_0*Sd_0*b_i^2*eps_s^2-Sn_0^2*h2^2+Sd_0*h1^2+Sd_1, 18335859437563935677633365-In_0, -Ad_0^2*f^2*g_ai^2-An_0^2*f^2*g_ai^2+In_0*dlt+In_0*g_ir+In_1, -Sd_1-2567645581378966339466914242565792044762080288781193540709483017994855271929822918049798678400525992150339367283492019695407974, (((Ad_0^2*eps_a^2+An_0^2)*b_a^2+b_i^2*In_0)*Sd_1+Sd_0*((Ad_1^2*eps_a^2+An_1^2)*b_a^2+In_1*b_i^2))*eps_s^2+Sd_1*h1^2-h2^2*Sn_1^2+Sd_2, -Ad_0^2*Sn_0^2*b_a^2*eps_a^2*eps_s^2-An_0^2*Sn_0^2*b_a^2*eps_s^2-In_0*Sd_0*b_i^2*eps_s^2+Ad_0^2*g_ai^2+Ad_0^2*h1^2-An_0^2*h2^2+Ad_1^2, -Ad_0^2*Sn_0^2*b_a^2*eps_a^2-An_0^2*Sn_0^2*b_a^2-In_0*Sn_0^2*b_i^2-Ad_0^2*h1^2+An_0^2*g_ai^2+An_0^2*h2^2+An_1^2, Ad_0^2*Sn_0^2*b_a^2*eps_a^2+An_0^2*Sn_0^2*b_a^2+In_0*Sn_0^2*b_i^2+Sn_0^2*h2^2-Sd_0*h1^2+Sn_1^2, -In_1+12260170483740471638265242403334153910409794874016812140454241890125755513060, (dlt+g_ir)*In_1-f^2*(Ad_1^2+An_1^2)*g_ai^2+In_2, -Sd_2+96752863147379284657867218516980119532068907690911332947363538119039273108099114588407531992691552742818522945876405893547769573341379853146485407670818505820728236634276696792634460420106446425791265540328587658661455722143434, (((Ad_0^2*eps_a^2+An_0^2)*Sd_2+(2*Ad_1^2*Sd_1+Ad_2^2*Sd_0)*eps_a^2+An_2^2*Sd_0+2*Sd_1*An_1^2)*b_a^2+b_i^2*(In_0*Sd_2+2*In_1*Sd_1+In_2*Sd_0))*eps_s^2+Sd_2*h1^2-h2^2*Sn_2^2+Sd_3, ((-Ad_0^2*Sn_1^2*eps_a^2-Ad_1^2*Sn_0^2*eps_a^2-An_0^2*Sn_1^2-An_1^2*Sn_0^2)*b_a^2-b_i^2*(In_0*Sd_1+In_1*Sd_0))*eps_s^2+(g_ai^2+h1^2)*Ad_1^2-An_1^2*h2^2+Ad_2^2, (-Sn_0^2*An_1^2-Sn_0^2*eps_a^2*Ad_1^2-Sn_1^2*(Ad_0^2*eps_a^2+An_0^2))*b_a^2+(g_ai^2+h2^2)*An_1^2-In_0*Sn_1^2*b_i^2-Sn_0^2*b_i^2*In_1-Ad_1^2*h1^2+An_2^2, ((Ad_0^2*eps_a^2+An_0^2)*b_a^2+b_i^2*In_0+h2^2)*Sn_1^2+Sn_0^2*(Ad_1^2*eps_a^2+An_1^2)*b_a^2+Sn_0^2*b_i^2*In_1-Sd_1*h1^2+Sn_2^2, -In_2+588858839417517666107400543095007538245504894115799254996153631410596688976478314349343541770790039499048342846952888805023272826639179537416450394449277856817683935512829623510, (dlt+g_ir)*In_2-f^2*(Ad_2^2+An_2^2)*g_ai^2+In_3, -Sd_3+13940907365117419012145526193663399631336478295503430577500876202612521225717790734668378536484768827804243062863497889704918204426947973934496828291973238998688399678158954741068733066557423711537108105318685851824281990625390540328458101084974988949433468131684842328375958597593549673302000365754683905838864724661039227768786, (((Ad_0^2*Sd_3+3*Ad_1^2*Sd_2+3*Ad_2^2*Sd_1+Ad_3^2*Sd_0)*eps_a^2+3*An_1^2*Sd_2+3*An_2^2*Sd_1+An_3^2*Sd_0+Sd_3*An_0^2)*b_a^2+b_i^2*(In_0*Sd_3+3*In_1*Sd_2+3*In_2*Sd_1+In_3*Sd_0))*eps_s^2+Sd_3*h1^2-h2^2*Sn_3^2+Sd_4, (((-Ad_0^2*Sn_2^2-2*Ad_1^2*Sn_1^2-Ad_2^2*Sn_0^2)*eps_a^2-2*Sn_1^2*An_1^2-An_2^2*Sn_0^2-Sn_2^2*An_0^2)*b_a^2-b_i^2*(In_0*Sd_2+2*In_1*Sd_1+In_2*Sd_0))*eps_s^2+(g_ai^2+h1^2)*Ad_2^2-An_2^2*h2^2+Ad_3^2, ((-Ad_0^2*Sn_2^2-2*Ad_1^2*Sn_1^2-Ad_2^2*Sn_0^2)*eps_a^2-2*Sn_1^2*An_1^2-An_2^2*Sn_0^2-Sn_2^2*An_0^2)*b_a^2+(g_ai^2+h2^2)*An_2^2-In_0*Sn_2^2*b_i^2-2*b_i^2*Sn_1^2*In_1-Sn_0^2*b_i^2*In_2-Ad_2^2*h1^2+An_3^2, ((Ad_0^2*eps_a^2+An_0^2)*Sn_2^2+(2*Ad_1^2*Sn_1^2+Ad_2^2*Sn_0^2)*eps_a^2+An_2^2*Sn_0^2+2*Sn_1^2*An_1^2)*b_a^2+(In_0*b_i^2+h2^2)*Sn_2^2+2*b_i^2*Sn_1^2*In_1+Sn_0^2*b_i^2*In_2-Sd_2*h1^2+Sn_3^2, -In_3+38886779554524151594439548908964211387191427476971839694922533500156000338220625138366341714147960155866436459516359636665355288577170331726361195146347616821790867107747682159281443996776543112641397595491033382475944921375469051285888128137688703995134150576597492136421333485, (dlt+g_ir)*In_3-f^2*(Ad_3^2+An_3^2)*g_ai^2+In_4, -Sd_4-689321116516347935982262570171484850683734204148402173860172133151410396697432705194647587827023645553733299012760428316566776030782871179238865061643553884952165029050354303602134844214527894243905044616748453097048719918810720467689621509420304847936106888552568102883624146950245808181686567043027949578578198152272947660787657643647247124715303250454639305915791020313066673170204247903595764523309453772613948664010517404586, (((Ad_0^2*Sd_4+4*Ad_1^2*Sd_3+6*Ad_2^2*Sd_2+4*Ad_3^2*Sd_1+Ad_4^2*Sd_0)*eps_a^2+An_0^2*Sd_4+4*An_1^2*Sd_3+6*An_2^2*Sd_2+4*An_3^2*Sd_1+An_4^2*Sd_0)*b_a^2+b_i^2*(In_0*Sd_4+4*In_1*Sd_3+6*In_2*Sd_2+4*In_3*Sd_1+In_4*Sd_0))*eps_s^2+h1^2*Sd_4-h2^2*Sn_4^2+Sd_5, (((-Ad_0^2*Sn_3^2-3*Ad_1^2*Sn_2^2-3*Ad_2^2*Sn_1^2-Ad_3^2*Sn_0^2)*eps_a^2-Sn_3^2*An_0^2-3*An_1^2*Sn_2^2-3*An_2^2*Sn_1^2-An_3^2*Sn_0^2)*b_a^2-b_i^2*(In_0*Sd_3+3*In_1*Sd_2+3*In_2*Sd_1+In_3*Sd_0))*eps_s^2+(g_ai^2+h1^2)*Ad_3^2-An_3^2*h2^2+Ad_4^2, ((-Ad_0^2*Sn_3^2-3*Ad_1^2*Sn_2^2-3*Ad_2^2*Sn_1^2-Ad_3^2*Sn_0^2)*eps_a^2-Sn_3^2*An_0^2-3*An_1^2*Sn_2^2-3*An_2^2*Sn_1^2-An_3^2*Sn_0^2)*b_a^2+(-In_0*Sn_3^2-3*In_1*Sn_2^2-3*In_2*Sn_1^2-In_3*Sn_0^2)*b_i^2+(g_ai^2+h2^2)*An_3^2-Ad_3^2*h1^2+An_4^2, ((Ad_0^2*Sn_3^2+3*Ad_1^2*Sn_2^2+3*Ad_2^2*Sn_1^2+Ad_3^2*Sn_0^2)*eps_a^2+3*An_1^2*Sn_2^2+3*An_2^2*Sn_1^2+An_3^2*Sn_0^2+Sn_3^2*An_0^2)*b_a^2+(In_0*b_i^2+h2^2)*Sn_3^2+(3*In_1*Sn_2^2+3*In_2*Sn_1^2+In_3*Sn_0^2)*b_i^2-Sd_3*h1^2+Sn_4^2, -In_4+2567986625823541032588811111380083067027795413568243618946870672877728124300872421736406240110047618969316285414744056672065949983325230523934183711276436767413944880796234134219386211836185377287435100523880513769276484168933501241729403975629191152927519614074758412725168407961070815921661303090522170708844072843202902759656028767796580638336954802550352964925207489549663985, (dlt+g_ir)*In_4-f^2*(Ad_4^2+An_4^2)*g_ai^2+In_5, -Sd_5-315077053352381591929029121295133280749299656650421770544899007307495523341635634035848072048819568659763828348498937513151968749978543420488587789125209367552210717766961346033881153353601802045908090652085503664256158611956419789715685118171580918828663677065897544508691889839391425639650338602368883588667712613123743402170714692098121506064475637948147070216786311284778757677255149907918959608082364038573668692714966822934149375195559275159685792058798233035499792978730678027272951017888870259967804827137783180582356354384, (((Ad_0^2*Sd_5+5*Ad_1^2*Sd_4+10*Ad_2^2*Sd_3+10*Ad_3^2*Sd_2+5*Ad_4^2*Sd_1+Ad_5^2*Sd_0)*eps_a^2+An_0^2*Sd_5+5*An_1^2*Sd_4+10*An_2^2*Sd_3+10*An_3^2*Sd_2+5*An_4^2*Sd_1+An_5^2*Sd_0)*b_a^2+b_i^2*(In_0*Sd_5+5*In_1*Sd_4+10*In_2*Sd_3+10*In_3*Sd_2+5*In_4*Sd_1+In_5*Sd_0))*eps_s^2+Sd_5*h1^2-h2^2*Sn_5^2+Sd_6, (((-Ad_0^2*Sn_4^2-4*Ad_1^2*Sn_3^2-6*Ad_2^2*Sn_2^2-4*Ad_3^2*Sn_1^2-Ad_4^2*Sn_0^2)*eps_a^2-An_0^2*Sn_4^2-4*An_1^2*Sn_3^2-6*An_2^2*Sn_2^2-4*An_3^2*Sn_1^2-An_4^2*Sn_0^2)*b_a^2-b_i^2*(In_0*Sd_4+4*In_1*Sd_3+6*In_2*Sd_2+4*In_3*Sd_1+In_4*Sd_0))*eps_s^2+(g_ai^2+h1^2)*Ad_4^2-An_4^2*h2^2+Ad_5^2, ((-Ad_0^2*Sn_4^2-4*Ad_1^2*Sn_3^2-6*Ad_2^2*Sn_2^2-4*Ad_3^2*Sn_1^2-Ad_4^2*Sn_0^2)*eps_a^2-An_0^2*Sn_4^2-4*An_1^2*Sn_3^2-6*An_2^2*Sn_2^2-4*An_3^2*Sn_1^2-An_4^2*Sn_0^2)*b_a^2+(-In_0*Sn_4^2-4*In_1*Sn_3^2-6*In_2*Sn_2^2-4*In_3*Sn_1^2-In_4*Sn_0^2)*b_i^2+(g_ai^2+h2^2)*An_4^2-Ad_4^2*h1^2+An_5^2, ((Ad_0^2*Sn_4^2+4*Ad_1^2*Sn_3^2+6*Ad_2^2*Sn_2^2+4*Ad_3^2*Sn_1^2+Ad_4^2*Sn_0^2)*eps_a^2+An_0^2*Sn_4^2+4*An_1^2*Sn_3^2+6*An_2^2*Sn_2^2+4*An_3^2*Sn_1^2+An_4^2*Sn_0^2)*b_a^2+(In_0*Sn_4^2+4*In_1*Sn_3^2+6*In_2*Sn_2^2+4*In_3*Sn_1^2+In_4*Sn_0^2)*b_i^2-h1^2*Sd_4+h2^2*Sn_4^2+Sn_5^2, -In_5+169583477622830147132968249577916550812144602246633473543761712308157581568655150596506528613687730174248214564996511382055236831499985460522860052356183208293944323947786682376024041545808564590838087554610920566652053291674661480312368925689035656200789969449600797590651652400958158803879624170176033052408661648671189076671456732061819877825903458271558096024446905814291906889995840215719096706632079645922598372593823060196059991697176874485212941661009852949465369819568185, (dlt+g_ir)*In_5-f^2*(Ad_5^2+An_5^2)*g_ai^2+In_6, -Sd_6-13571001718718277371759931425380767555419156751413550411308400675287883192151324622273248550769977512828383454793986726536723067964655756716152276774859488278637919422767064454507771189824109744922803910491716818354592874516354105166863063518296190418780878425379036004832028109865093397620638664864833954752055866050483896417200411530148380687652112874728134228067688787702429557737702943950883468939699983889150472317693665353543640360337437857437486715372772623176769156258971148836470306615746822421894147066997294538733389013850422311828295367451032575122082720091534155048933864353395412132584378217773401529719500320206076396, (((Ad_0^2*Sd_6+6*Ad_1^2*Sd_5+15*Ad_2^2*Sd_4+20*Ad_3^2*Sd_3+15*Ad_4^2*Sd_2+6*Ad_5^2*Sd_1+Ad_6^2*Sd_0)*eps_a^2+An_0^2*Sd_6+6*Sd_5*An_1^2+15*Sd_4*An_2^2+20*Sd_3*An_3^2+15*Sd_2*An_4^2+6*Sd_1*An_5^2+Sd_0*An_6^2)*b_a^2+b_i^2*(In_0*Sd_6+6*In_1*Sd_5+15*In_2*Sd_4+20*In_3*Sd_3+15*In_4*Sd_2+6*In_5*Sd_1+In_6*Sd_0))*eps_s^2+Sd_6*h1^2-h2^2*Sn_6^2+Sd_7, (((-Ad_0^2*Sn_5^2-5*Ad_1^2*Sn_4^2-10*Ad_2^2*Sn_3^2-10*Ad_3^2*Sn_2^2-5*Ad_4^2*Sn_1^2-Ad_5^2*Sn_0^2)*eps_a^2-An_0^2*Sn_5^2-5*Sn_4^2*An_1^2-10*Sn_3^2*An_2^2-10*Sn_2^2*An_3^2-5*An_4^2*Sn_1^2-An_5^2*Sn_0^2)*b_a^2-b_i^2*(In_0*Sd_5+5*In_1*Sd_4+10*In_2*Sd_3+10*In_3*Sd_2+5*In_4*Sd_1+In_5*Sd_0))*eps_s^2+(g_ai^2+h1^2)*Ad_5^2-An_5^2*h2^2+Ad_6^2, ((-Ad_0^2*Sn_5^2-5*Ad_1^2*Sn_4^2-10*Ad_2^2*Sn_3^2-10*Ad_3^2*Sn_2^2-5*Ad_4^2*Sn_1^2-Ad_5^2*Sn_0^2)*eps_a^2-An_0^2*Sn_5^2-5*Sn_4^2*An_1^2-10*Sn_3^2*An_2^2-10*Sn_2^2*An_3^2-5*An_4^2*Sn_1^2-An_5^2*Sn_0^2)*b_a^2+(-In_0*Sn_5^2-5*In_1*Sn_4^2-10*In_2*Sn_3^2-10*In_3*Sn_2^2-5*In_4*Sn_1^2-In_5*Sn_0^2)*b_i^2+(g_ai^2+h2^2)*An_5^2-Ad_5^2*h1^2+An_6^2, ((Ad_0^2*Sn_5^2+5*Ad_1^2*Sn_4^2+10*Ad_2^2*Sn_3^2+10*Ad_3^2*Sn_2^2+5*Ad_4^2*Sn_1^2+Ad_5^2*Sn_0^2)*eps_a^2+An_0^2*Sn_5^2+5*Sn_4^2*An_1^2+10*Sn_3^2*An_2^2+10*Sn_2^2*An_3^2+5*An_4^2*Sn_1^2+An_5^2*Sn_0^2)*b_a^2+(In_0*Sn_5^2+5*In_1*Sn_4^2+10*In_2*Sn_3^2+10*In_3*Sn_2^2+5*In_4*Sn_1^2+In_5*Sn_0^2)*b_i^2-Sd_5*h1^2+h2^2*Sn_5^2+Sn_6^2, -In_6+11198872919920368308801430700525371860907479083685623424597756190981661506872391307888896562847788232520383062456846067244163117967406529112552982157936041358366579142215899595714792293730381037809521032506876515887220925953515370051690983493527304514355522426632249791248041185031235969328922764300625529545003496087885747296226002109994685386936893081090175613967286870561849801571301056739133997899730961349397812043793865787073196062616213268538344619924159667263773505546705255008449901477005861126880012986386825851211280395721322951312594824128704975741439237235645813977760, (dlt+g_ir)*In_6-f^2*(Ad_6^2+An_6^2)*g_ai^2+In_7, -Sd_7+9635758537313794817296028972249786827371296205399345868922484040760492044268797858878970434620369144772224637555315222639942161081570373051037265039555282501497809175117774280734721227251941200961959188718357866917196872843725345955284896438964843532037768017645691599124943465181420644907847956080774808820252264754142002348616076712870342516907402991269980969082330030092889626707634456774698462975296818283159980440885915590158357266933946721588780877609325197925009807810478395938325695091281756105793052868523239418469456794520335231563208878040706429602070973662621539595151673609185506154153123423428976914934726463190034227180609716784332893605974482292599361525004163457467569294579547824173287321722675541198693623011781736, -In_7+739545835682531231329843844086739594837806935259158330446155483991162296094383311045197781223636456279201153609449383436645694377141336949915977344051328937684876581269100158543803559205318208750769432039195388197850301345390264307752002093734681007200731263238654155152344263734544288145642101404248312006961167700610933041979617034873436089167751885932363419624340257964404953173590258397466890275658836957786560437811580713335121448831364233876250632545019730931935245191000463062466867990654020649363920427654931499657465355445587628912324892355063803902089532560480915453579997181323412185341362261308578760345201613349840104147552985658380103005804954073343564091818162439060, z_aux-1];
vars:=[Sd_7, In_7, Sn_6, An_6, Ad_6, Sd_6, In_6, Sn_5, An_5, Ad_5, Sd_5, In_5, Sn_4, An_4, Ad_4, Sd_4, In_4, Sn_3, An_3, Ad_3, Sd_3, In_3, Sn_2, An_2, Ad_2, Sd_2, In_2, Sn_1, An_1, Ad_1, Sd_1, In_1, Sn_0, An_0, Ad_0, Sd_0, In_0, z_aux, w_aux, b_a, b_i, dlt, eps_a, eps_s, f, g_ai, g_ir, h1, h2];
new_weights:={Ad_0 = Ad_0, Ad_1 = Ad_1, Ad_2 = Ad_2, Ad_3 = Ad_3, Ad_4 = Ad_4, Ad_5 = Ad_5, Ad_6 = Ad_6, An_0 = An_0, An_1 = An_1, An_2 = An_2, An_3 = An_3, An_4 = An_4, An_5 = An_5, An_6 = An_6, In_0 = In_0^2, In_1 = In_1^2, In_2 = In_2^2, In_3 = In_3^2, In_4 = In_4^2, In_5 = In_5^2, In_6 = In_6^2, In_7 = In_7^2, Sd_0 = Sd_0^2, Sd_1 = Sd_1^2, Sd_2 = Sd_2^2, Sd_3 = Sd_3^2, Sd_4 = Sd_4^2, Sd_5 = Sd_5^2, Sd_6 = Sd_6^2, Sd_7 = Sd_7^2, Sn_0 = Sn_0, Sn_1 = Sn_1, Sn_2 = Sn_2, Sn_3 = Sn_3, Sn_4 = Sn_4, Sn_5 = Sn_5, Sn_6 = Sn_6, b_a = b_a, b_i = b_i, eps_a = eps_a, eps_s = eps_s, f = f, g_ai = g_ai, h1 = h1, h2 = h2, z_aux = z_aux^2};
gb:=CodeTools[Usage](Groebner[Basis](subs(new_weights, et_hat), tdeg(op(vars)), characteristic=11863279),output='all');
# {Ad_0 = Ad_0, Ad_1 = Ad_1, Ad_2 = Ad_2, Ad_3 = Ad_3, Ad_4 = Ad_4, Ad_5 = Ad_5, Ad_6 = Ad_6, An_0 = An_0, An_1 = An_1, An_2 = An_2, An_3 = An_3, An_4 = An_4, An_5 = An_5, An_6 = An_6, In_0 = In_0^2, In_1 = In_1^2, In_2 = In_2^2, In_3 = In_3^2, In_4 = In_4^2, In_5 = In_5^2, In_6 = In_6^2, In_7 = In_7^2, Sd_0 = Sd_0^2, Sd_1 = Sd_1^2, Sd_2 = Sd_2^2, Sd_3 = Sd_3^2, Sd_4 = Sd_4^2, Sd_5 = Sd_5^2, Sd_6 = Sd_6^2, Sd_7 = Sd_7^2, Sn_0 = Sn_0, Sn_1 = Sn_1, Sn_2 = Sn_2, Sn_3 = Sn_3, Sn_4 = Sn_4, Sn_5 = Sn_5, Sn_6 = Sn_6, b_a = b_a, b_i = b_i, eps_a = eps_a, eps_s = eps_s, f = f, g_ai = g_ai, h1 = h1, h2 = h2, z_aux = z_aux^2}
quit;