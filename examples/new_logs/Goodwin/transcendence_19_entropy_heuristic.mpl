infolevel[Groebner]:=10;
Et_hat := [2770176839777965596-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 4274090482076395323-gama_0, gama_1, -x1_1-82435295849868643293560130417543511343624170707651846599/11687513682358195325, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+381416772179959954766741358753654394130874372014387140569425596766510290644999644445287839549884633408229132824/21238563697205592365409768497229898485248548024549727275, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-1764762928310777543321898573287689873954182096764771640565815607971133214134264399693884536354779547455164980799646965800655392087391844676168049532977612051606997674/38594743089040417549580624943330815625253922837275467017917778726118488895903945737712610925, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+122516625855618242496385513652441182491026142911482932922977050062580917197572045779854564293707588290598512633100643761827636229639928918710270035288633774826997768568672123970893934407339265206227464149727246625024714019/70134412823076977208492363760722134014026897395599179088419388745233242982684331632794229208422675773333749791832877970415916475, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-3246920517540898058049595017427122488445290269049296824859161197461285966715144318216431001411272086914152220545676294002659611533504657905387747855810984623427863988734772376567323726719395206500956754963566670195906859821325510825054912806359948088726484747908607833832180854775498130024429329/127448337994885238668413320616433161720768811172469359935670460069766129632219511887841053236507210023427056703363684025449421854587234585521942383944604544336895325, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6-1251119301165192605861418021733696118638443282715980547907239346975757228891240868330242638439748741392445973799633940172262581667699946552934772090044623988942467563200788655993322486397231332517375383265282312213827057718719609928966508917625376941114624783325513937575068440211579473251441005427425981900481212108992706438823611938209912249951032357770971140883559281/231599270655245255687901522199527455407477906580065449781246422720034402288307036201650581861097411714343663288526101036244584975802101512634913253652646062411040774597724402458670581661154576960627275, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+3087863808451789393753571450864037993132641729832290648640004252517512944917203121707977028623001492091346184801973690464369326806551696294974105870473069533588747946866121976052025483306496603442006631735244555402681050466085109153145334742140183420458387795292061035071509786615302288833662137109167432525112560432612528212274655315676367252451467491061756019011922366460793913001232437048946971305885522297684598300609280857890391391278264/16834498750448721359409178220756050263203971911627632941696909315752397742808958862590782569986734857764480015509008244909421639676315793987645283475782741586474326876424424451147497531675025640558103602842739695029799006750454500356437, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8+4090642138277790995344918494251611647965979181187087831171909236253071375585234258951127113107280957843991670229581910429866384125169824627443301219720097307755772500549239571349872109850754123482922548741778110591271048003268679287986759221041064965391361010569209580732494393456860954893526343834384978734673362018327342089087702375646230260812717061255138664078625427688217424811426689332251103448178025095731811418133747302568039877825276263224226493793768803812422687935717211299447746928271763182317439916701149/109255968053699603295867040844574198621055350619295482560467167632453559852655020981828759210241869410503720590798231649156748459213674955036279423949653050184299286957912177287758277500648801451613420779442888256866964046033870084339536700976068267351040327082135061430925, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9-2239170044441205669307143460530324778646484747012836355748557226861469459036338263050757903449125386426892487087787653706114588898653919260833399304855525236927537144202488549983488475932522726733615110108500227674048656633168196793677342760258360424736149748561976287870340791299459926702367085195789240773524519624352528544200292531162793858662903057923539295799807999334692705697813038304874769819850791167690523267802105916245409337540481324190085990944156021203331973626055380506394368781906286822263938191794912804021085920574805112552042598975847691164702351077558321043556790646722738/992700372325965080103002472165611611731553766534029729503243124201657607329506459017081491210887106345097400137181880640051298882799164118478834887645566892532539737261128135774314524462627620180765642880253004521940196826554355890250863722254623402760033544893551280811537607309868601213295081197465368282375, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -9218741073459639946587347981279250174-x2_1, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x2_1] 161
# 274, [1 = 7, -1 = 3, -sgm*x2_2*x4_0 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -6*sgm*x2_2*x4_2 = 1, x3_6*x4_1 = 1, -7*sgm*x2_6*x4_1 = 1, b*x1_1*x4_2 = 1, x3_1*x4_4 = 1, -sgm*x2_7*x4_0 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, -6*sgm*x2_5*x4_1 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, -sgm*x2_0*x4_4 = 1, x1_5*x4_2 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, -sgm*x2_6*x4_0 = 1, delta*sgm*x3_1*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -sgm*x2_0*x4_6 = 1, x3_2*x4_3 = 1, -x2_2 = 1, -sgm*x2_0*x4_7 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, -20*sgm*x2_3*x4_3 = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_0*x4_2 = 1, -35*sgm*x2_4*x4_3 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, -15*sgm*x2_2*x4_4 = 1, x1_7*x4_0 = 1, -5*sgm*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, -3*sgm*x4_2 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, -15*sgm*x2_4*x4_2 = 1, b*x1_5*x4_3 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, -21*sgm*x2_5*x4_2 = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, -x2_4 = 1, x1_3*x4_5 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, x3_4 = 1, b*x1_1*x4_5 = 1, -4*sgm*x2_3*x4_1 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, x1_8*c = 1, -sgm*x2_4*x4_0 = 1, x1_2*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, -x2_0*x4_0*sgm = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, -sgm*x2_0*x4_3 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, x1_3*c = 1, -4*sgm*x4_3 = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, -sgm*x2_0*x4_2 = 1, x1_1*x4_8 = 1, -x2_5 = 1, -sgm*x2_5*x4_0 = 1, x2_4 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -x1_1 = 1, -x2_0 = 1, -x1_6 = 1, -sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, -x2_3 = 1, x3_1*x4_6 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, -2*sgm*x4_1 = 1, x1_3*x4_6 = 1, -sgm*x2_0*x4_1 = 1, x3_0*x4_4 = 1, -x4_0*sgm = 1, x3_1*x4_5 = 1, x3_3*x4_5 = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, -7*sgm*x4_6 = 1, -5*sgm*x2_4*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -x1_4 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, -35*sgm*x2_3*x4_4 = 1, b*x1_7*x4_1 = 1, -10*sgm*x2_3*x4_2 = 1, b*c*x1_4 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, -3*sgm*x2_2*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, -21*sgm*x2_2*x4_5 = 1, beta = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, -sgm*x2_0*x4_5 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, -x2_6 = 1, x3_7*x4_1 = 1, x3_1*x4_3 = 1, x3_3*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, -10*sgm*x2_2*x4_3 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, x2_7 = 1, x3_3*x4_2 = 1, -6*sgm*x4_5 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.581916971
# 2
# [x2_1 = 2, gama = 43]
# [x2_1 = [1, 2], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]