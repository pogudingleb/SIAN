infolevel[Groebner]:=10;
Et_hat := [3850727958025885559-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-77995420055065956401813935766921413921252294159937436839/10720774720679606970, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+102112049734941578430993478132894522104230149702951802736619179102932793833266054175733070711797508346546032437/7429061775949728092993200258534262850897787136332891700, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-401057037466160935665740368637665416519560456789209648254026086455831122064196975349222581147191298718636382661338082936903292407139500808949166157677640379979860851/15444114900880601264922691225808285904011211924831960345947878663267889044216999652163011000, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+966934662217644024544429410655495952332679389864607472984247902565184267900711854123155190946981767181110488741072951221438785600568017609762109971324294569349737726352444814383692037880593485709608435122526988846403819801/96319303405352222200364735504253097347797354971993681004583028927396342419230101059090438728485402006381721085590174936418390000, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-655542624035097504294475678864257153973987704497642301094568001987995602966705795684458389333467253725164567360230526085240926244916604212572990966113752518536855985006786755420106879989306730083372435848153051696108416647910659662628909050836753193633244263824797613206580330493303942976112644913/600708313039247829455024436002682627957222609242062776091708112036819941359605249554170816664399476663341735357294064875704943256940180069133686053862354610081100000, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-15491068429932373586436322776399304570424961785713658005048322078781141488137839848394861196292502369647898253806416003597262044220681366206049851579349059835424351706249784177696822019113853001344741821217554511476443692386843642127867663453953836661856257785615642533412412673371065198512773494682501229155461056949701189732258988168968605857371860643682000783150884811/249759888618136474558349061493629674217129131758082769159989979489463394026272367390445858719665487785358267459215558935423216626178421346304459332379172428550075297811965947647005346373924515602600000, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+34585469358650939551044349812646989237836536180681208852603216472622162078928025506109546122392594421384531553270272727833987360428579838851980881294942286929325152909054204797250901401813626566233296741351970983650982976439408874276943770188879943516716782088806547707076502501816221281484537955628585667548060020179390842908641652742854481842262330342269790361265260259697217607993122624401780224889383350441733809526003649889208073829674360371/1557661196170302115309052615070678851104920188645375897354500420068719698765663290846937920407407182176213181904190671816916363663093061552691677607949951531460610686523221608662675606513533785116126689736782105496392210679296319074000000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+160284209749499415282693570745614480099977409372691120406194775388559017849589422794198062595158521235728079720902396583372589916684426248029004824489230006264294668729711464196292427941575289895100044839257959068842538306567383019391474901735266871395337558687404542856406878782155168277621091274946075949812611313877930010035469490535267986262717974869142446715999106396581244515110113968829705657383339723254797903845570564475769526049687523043661611571356117327628130936483675463159960587706574492234251523623194257/3238187972561564299050428339180861286808362083496747500565829639882830707682153887197558722966283041165423388862273885167022129772269514701773413107898762463901311961956476028478607338054857399892208169839578404457842884603338798658375986630454774189867376776291969420000000, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-14486213391930939585033463098167724307272240517869058935158318177771866885904905071708102957139671286568891066204856600476740829250192128128852058120547173134709570231644563282695592274447617891169070419159130444129827264384362829386123439733374380474877157228152397970514071969661347385346711076299494600086579521819255111346972704256276518994157470148933100617803757811151635909228325600689418394386514846460099874482719660337688110152741475194188804802995220836444798665346758251849411045426163030153069313736530200652365099554939948552009216215404494261392522525636155513000936392881209896817/20195395580418506052477889305918069397191048662234461486171017800071471205358228239376426519533493421445262717867182530926694194382386960077498297998397958017846754335682046103740712842055844520125756809570802411214330680039829967829863217734402827363401619058483496350567270020771492779070793479095475800000000, 7177678245210539294-x2_0, -11124114685380807000331426998713649331203361107071532084704781141938171858969843269288076217569382797965825909269224336390570184964443517227923584415943327996387544058144193089797467870758016872326328430459314957227995408542422634188213679725858253877064003507447220698896215580285083195511591144788421446251890820071172936290930024377/600708313039247829455024436002682627957222609242062776091708112036819941359605249554170816664399476663341735357294064875704943256940180069133686053862354610081100000-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_0, x3_7] 43
# 275, [1 = 7, -1 = 2, -7*gama*sgm*x2_6*x4_1 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_6*x4_1 = 1, b*x1_1*x4_2 = 1, x3_1*x4_4 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, -gama*x4_4*sgm = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*sgm*x2_5*x4_0 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -21*gama*sgm*x2_5*x4_2 = 1, -15*gama*sgm*x2_2*x4_4 = 1, x3_2*x4_3 = 1, -gama*x2_4 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, delta*x4_0*sgm = 1, x1_7*x4_2 = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_7*x4_0 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, -gama*x4_5*sgm = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_3 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, -gama*x4_6*sgm = 1, x3_4 = 1, b*x1_1*x4_5 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, -gama*x4_1*sgm = 1, x1_8*c = 1, x1_2*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -5*gama*sgm*x2_1*x4_4 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, x1_3*c = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, x1_1*x4_8 = 1, x2_4 = 1, x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -gama = 1, -gama*x2_1 = 1, -x1_1 = 1, -x1_6 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, -gama*x4_7*sgm = 1, x3_1*x4_6 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, x3_0*x4_4 = 1, x3_1*x4_5 = 1, x3_3*x4_5 = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -x1_4 = 1, x2_1 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -gama*sgm*x4_0 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*c*x1_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -gama*x2_5 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, -gama*x4_3*sgm = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, -gama*x2_2 = 1, -gama*x2_6 = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -35*gama*sgm*x2_3*x4_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_2*x4_1 = 1, beta = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_1*x4_3 = 1, x3_3*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, -6*gama*sgm*x2_1*x4_5 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, -gama*sgm*x2_4*x4_0 = 1, -gama*x4_2*sgm = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, -10*gama*sgm*x2_2*x4_3 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.588688399
# 2
# [x2_0 = 10, x3_7 = 3]
# [x2_0 = [4, 4, 2, 2, 4, 4, 4, 4, 4, 4], x3_7 = [4, 2, 1]]