infolevel[Groebner]:=10;
Et_hat := [2101659478933730186-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-43741664279997391663504342265054944406188785384621108429/11780846882724056455, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+18954424021620502033449835638524157731675588931861748769706680104089154615621380556864191366357380332203152613/2889584237224052773741145753155494661286031996110264905, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-3733388623919897652368767279998989183210979578316407291394370618182987002496948749814819403235449922378866570115116954001939203034166010453184897882061949767385197/322159933412306403886557071894605342497603360891288891141690371399494256087931437804754025, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+261144057426930646439349891144428752851644955804000127200409180415232346507092819004142886934790467143020725073522055357165118245752725038025100602488142675424680714147245923963145138207322970959608271155852676606627619937/4346033448156951714440277863569875903114095882785901668012504992470792509877398579738319047629747612458701078341704001152082625, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+2756797774810165536538893451173002831830323057064785460319618060965502549937409183580193675625255552896673227264719749402846706212659242102314569851603797457253100762577966659385907495841572893533769254303043864480431561328019291208525076152412884973709079713127794193263347259969146555491728057/5329935051044062304603455942055674065682031716491159534364047481138852115983364296248275847298252963051697787971251420469239296401969321803167662832111928891656875, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-15864852469999834791484589736933699874555442721894256092395760310755723191497163881928410746736189599922936054568752970690175273800448390807619588932750453514393126998642377937951879912555597745215161987341984680748360482488211484772004070543505880344241212835244910683313673166743329254337958486579349020645669326844797846725489471982045630217498227342648064226698009579/6536582837482603553479542306929118060925481355176348995318935435211106604167481438028597696565495550914492968512816519241337285189424060798386283095645128477831937853499728165703316725615985408565625, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-71243251500227521604662955314842362426693372223935712821131714526302889279470843505553670685617944938218991398065631471040703797796242578358407399371667002779187484430507462499670122435230728417057372504334784193646674957595193125215199794800059886029449846979287916430317630113148091332427924996820551302730521199468139324153753096409397511369524297230310904329881187767592300894411397243151147754431249105066039492093429405864905337157304383/1603280894873287397698837019845163806867583998831753040611525634170427438587559305984515439265563319973300968884816524303048187531726222692967284167771833913150777827199816130318500039190995027601311268989797372389028662682468899784375, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+37237255463486940139201764217831596194059084568916874728143459198913477360551028939720375111666880957819555927703199920906444885405213595263007928277051018874927315336537361493897717597021518973059105350990596202136086171088203627657026632946545949608032625111321175000475481770213047406726183332545753800169375938958900057614538501665803862513455003000271717050567443944625941223695070267159107721990957094561252354202567666606444765112672111437222256834370000814937693268990782410992107292526165255149999794376118711/178749890617413373894616141273136442925054453753451422940098835189639733579535428090209052832149556594126994787076261277761859175077399156087745691257999973436853416508640560876016285339617611958715297923558746432422704948554110518186856883417429711847516843181970484375, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+15797299106882229457764242839702374329538874417475356985708538348449834704808655927236205131563273729214081952116294713522356194506212329962456430622374186211915509825696435122779406168141649733821951218478013473297696360571866091962033323589224608280313776104484911877236389233123243042572030602171579778885896965272153000237126920503965185242838556016357051263185811375745230813657258437451152444839310283210555028539437057934318536091362485168946415587689134686891160837737765832936393373692181976432454519702990606080721948070099197627471982747817572789341176030794916862931920291049900209/2411389260139446831081213756420915061117454132398708146040950545308800853603494464481478193676698811112511126770366238487715498368795236220980962853742823444661401642436709028790766319996506034194243800418624703750627993595722341207860723035000428419347359374956741302618774875404673740457348859011532109375, 7153050780407862647-x2_0, 18483960468859377656024139699934615081-x3_1, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_0, x3_1] 81
# 275, [1 = 7, -1 = 2, -7*gama*sgm*x2_6*x4_1 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_6*x4_1 = 1, b*x1_1*x4_2 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, -gama*x4_4*sgm = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, x1_5*x4_2 = 1, delta*sgm*x4_1 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*sgm*x2_5*x4_0 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x4_5 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_6*x4_3 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -21*gama*sgm*x2_5*x4_2 = 1, -15*gama*sgm*x2_2*x4_4 = 1, x3_2*x4_3 = 1, x4_3 = 1, -gama*x2_4 = 1, delta*sgm*x3_2*x4_3 = 1, -alpha*x1_3 = 1, x4_2 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, delta*x4_0*sgm = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_7*x4_0 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, x3_2 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_3*x4_4 = 1, delta*sgm*x4_6 = 1, b*x1_0*x4_1 = 1, -gama*x4_5*sgm = 1, beta*x2_4 = 1, delta*sgm*x4_4 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_3 = 1, x1_3*x4_5 = 1, x4_4 = 1, -3*gama*sgm*x2_1*x4_2 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, -gama*x4_6*sgm = 1, x3_4 = 1, delta = 1, b*x1_1*x4_5 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, -gama*x4_1*sgm = 1, x1_8*c = 1, x1_2*x4_3 = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_2*x4_2 = 1, x3_2*x4_5 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x4_5 = 1, b*x1_6*x4_0 = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*c*x1_1 = 1, x1_5*c = 1, x1_3*x4_1 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -5*gama*sgm*x2_1*x4_4 = 1, b*x1_4*x4_3 = 1, delta*x3_2 = 1, x1_3*c = 1, x3_2*x4_1 = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, x1_1*x4_8 = 1, x2_4 = 1, x4_1 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x3_2*x4_4 = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -gama = 1, -gama*x2_1 = 1, -x1_1 = 1, -x1_6 = 1, delta*sgm*x3_0*x4_5 = 1, delta*sgm*x4_3 = 1, b*x1_3*x4_0 = 1, -gama*x4_7*sgm = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, x3_0*x4_4 = 1, x3_3*x4_5 = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -x1_4 = 1, x2_1 = 1, beta*x2_2 = 1, b*x1_4*x4_0 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -gama*sgm*x4_0 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*c*x1_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -gama*x2_5 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, -gama*x4_3*sgm = 1, b*c*x1_5 = 1, x3_2*x4_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, -gama*x2_2 = 1, -gama*x2_6 = 1, delta*sgm*x3_2*x4_5 = 1, b*c*x1_6 = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -35*gama*sgm*x2_3*x4_4 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x1_4*x4_5 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_2*x4_1 = 1, beta = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, x3_2*x4_6 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x3_7*x4_1 = 1, x4_6 = 1, x3_3*x4_3 = 1, x1_2*x4_2 = 1, -6*gama*sgm*x2_1*x4_5 = 1, -10*gama*sgm*x2_3*x4_2 = 1, delta*sgm*x4_2 = 1, x2_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, -gama*sgm*x2_4*x4_0 = 1, -gama*x4_2*sgm = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, -10*gama*sgm*x2_2*x4_3 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, x4_7 = 1, delta*sgm*x3_0*x4_6 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.588688399
# 2
# [x2_0 = 10, x3_1 = 16]
# [x2_0 = [4, 4, 2, 2, 4, 4, 4, 4, 4, 4], x3_1 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2]]