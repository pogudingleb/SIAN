infolevel[Groebner]:=10;
Et_hat := [2769674076132215873-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-1483784466875697798703866130381954361595001845520713767/934083667105027188, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+1617739164954054276812961110425719892248012180591282971004080776023648793409128904224317445318807403887/1775689725555046029244738123413054515429372316246, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-1763787170071057102734630648825455721203318384802854920772011317445569671857208174696071872257125224181864253718519274234150259855700998855849414149075/3375579846304310798151579425157789822689570771708793298691215346389552101278157, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+43612037985084251373435190731154405190878943653168794637698147167575919190476119153503142327063771745123326150298898711608752858194050195570887306969549399032083107096616479665852011808290678698653998873895443/1501569873126401561919561882877845461688506824784986024715745889110924265138665060919434333556160570835377964471, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-982295288656833322415478455800292789592744909539956316966277886669863551990923827784481159247195504240935860896836091820002476325397797490758379901879057637687120559591499373464138234282611962674088578555222150664106002958404920485135667217716209856864910341437413859018333736572/74216453552396213498311130851417259010389640433202727991761315388228350212738515353124999282317826608675784643209887051425291982955479450698157, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-35858297864746450588063783918293606772487175894181059985859938497747212533669633295293052291603080100313941566270305733514845451940205665966359759635197533450446332441946606328617752439979410500803175096318220000552441159165166569805948266788918575383005635306319791390018050636464708885856510481685915103518593496918495184674156971838340172970783691/3668215563240270445570606036777164844831271706556145889074590603877661096965440519551156725336918241374616563785031320535726007354992808832443027384181178771760308051691802719, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+225511192059351531460184342999031583371925513561446414343721340183680685186735994857357812114332369667590853811543425115478210919025311241640943382605986682391350531508812916216713407405149631098276708561665626931900182740303900051322946515800826905401672914938669546493735569169554383301948016971296207664717581996323968370682207583807737493352237235086035402432262072483527611818906680592441749119651368955413200319114/20144986985636615221596044584140492902670392454608541130715157172092715025001826323835645225492616052819456195518299076221537455456688540765345270216681207325683503755448068985006477641577562925950271606397, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+153093897978133028588577705677011345349667869781250003859569856487756425168091070548038006639806534666424806551600025002764596214494607231274946500894093667837705855646653923727757507722482038917500184521360197454093212362995535321046291205214282536084505420211277533408363109658392858332152258481959729059283991758573093686705987626236198684316710811687157529991216741375497272578145608691249808614848811721606273307970035784326815147203341333462952861818086117656361142896833948205010430737/26883475073421568571050668450955305688057092814870556751443677275545792219202020379334407948144571974604046024697563101906940063145706436316798654571940299821614115254154551492389354051577539660020287799415823532327231007601936436691647573, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-62617107138961066384989422604087820588871403186439026855331503712220551273751745173020061207325964292910132940401749784649552225390908320622777118476886073178540636061689868108365104944069756688563218092265883730282754708367049837045491491241404081152299045899561667416154968576454570420952426993319256613184361869380259040173141802753849945743261703450905033428159921130619929043901930228322842447714671000029321889652232181540079475045533248895224851208597919066698742271150779280060278247974271833447648208280791558954175254899284444617304902792418737320084344/3986220451318894566433197884236300314033924642804448801778712804736039166231865258399013727116288659664043232955855723786213107890541202080483565593840214782682356031753576539586530636512575147256447204124590186498191547929417243505300470579976273142277155382336638837573, 32360277980051467821134091906433767151625677544065736792322088097960216705/934083667105027188-x2_2, 11177034950835434042655232129497119450861232694395094821672002477713953356939524984069134752/233520916776256797-x3_3, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_2, x3_3] 62
# 274, [1 = 8, -1 = 2, -7*gama*sgm*x2_6*x4_1 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_6*x4_1 = 1, b*x1_1*x4_2 = 1, x3_1*x4_4 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, -10*gama*x4_3*sgm = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -gama*sgm*x2_0*x4_1 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*sgm*x2_5*x4_0 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_0*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x4_4*delta*sgm = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, x4_2*delta*sgm = 1, -gama*sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -21*gama*sgm*x2_5*x4_2 = 1, x3_2*x4_3 = 1, x4_3 = 1, -gama*x2_4 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x4_2 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, delta*x4_0*sgm = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, x1_7*x4_0 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -gama*sgm*x2_7*x4_0 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_3 = 1, x4_4 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, x3_4 = 1, delta = 1, b*x1_1*x4_5 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, -6*gama*x4_2*sgm = 1, x1_8*c = 1, x1_2*x4_3 = 1, x4_3*delta*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, delta*sgm*x3_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x4_5 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -5*gama*sgm*x2_1*x4_4 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, x1_3*c = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, x1_1*x4_8 = 1, x2_4 = 1, x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x4_1*delta*sgm = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, -gama = 1, -gama*x2_1 = 1, -x1_1 = 1, -x1_6 = 1, -21*gama*x4_5*sgm = 1, -gama*sgm*x2_0*x4_7 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, x3_1*x4_6 = 1, -gama*sgm*x2_0*x4_5 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, x3_0*x4_4 = 1, x3_1*x4_5 = 1, -3*gama*x4_1*sgm = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, -x1_4 = 1, x2_1 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, b*x1_2*x4_4 = 1, -gama*sgm*x4_0 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*c*x1_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -gama*x2_5 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, -gama*x2_6 = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -35*gama*sgm*x2_3*x4_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, beta = 1, -gama*x2_0 = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x3_7*x4_1 = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, -6*gama*sgm*x2_1*x4_5 = 1, -10*gama*sgm*x2_3*x4_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, -gama*sgm*x2_0*x4_4 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, -15*gama*x4_4*sgm = 1, -gama*sgm*x2_4*x4_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, x2_7 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.577999884
# 2
# [x2_2 = 9, x3_3 = 12]
# [x2_2 = [4, 1, 4, 2, 2, 4, 4, 4, 4], x3_3 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2]]