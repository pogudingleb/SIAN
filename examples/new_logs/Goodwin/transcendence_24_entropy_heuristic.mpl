infolevel[Groebner]:=10;
Et_hat := [4217285159882646357-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 7380773711794010304-gama_0, gama_1, -x1_1-248586588054339700373073627066911185255521009790104726257/14584515611471787902, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+1774834447873695779455181086899390929420462766986036569688182638166248104413699153929439955951160519212680382279/25764365370022202012431836897089557793357510485106452502, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-12671791113165547002369779799518445802756472687179193775286797444111773989217269780178918194711729017801151698813908246420981668386015688683259069996591576381973254673/45514197427158297873559315279437244932403837664660565742524392698740329080394380397746697102, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+837991647502814081484385197998959276021509913582662188607600437332693095142289394521269086335943564779856539846259120810107847403308040033721276004378369298508078030992416669021290454567592598063453442740139142021282185431/80403384197021962382317625700673137307375970478231149277412329130131695930433508092528598487293806826761191741821718014326521702, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-26076759480627725415888620304773237502446288660024649976767427443935616981789518556048732930545646114714287838091810937759350836445398140805352935938051044542485736670099134701671099468742256826937926886287980330453487547946289717908556267314196819637189019390318526096133865109142143976514180737/142037090749104045826320283521355186997827349718655503748712482139921321719476130783244216116971319781664630958835510194424427159364423151859067573511374677979926302, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6-29542077533795369613170450101468417650724762406373674708986003863195201143087775276452128329394126005250906033161336510656887456802925253558821034674664940887602774580988485450421583653791324381368846703197618880420797625205195037871174565054574918121066528534551600606127622926411707245944093277358245415519909267668293595476943703388748007700781308177734426586919794201/250916492507742688616037181532939106391557309294198078292839406152347001775038431663421767278094671308941477041924615619791908771209753798688784594083070053325063345338436810951830390156229270040910902, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+2210703584804029050927766628135469877761029304897727010033859949821934234969785657659420862629441147104113666716537596352956925461554263880801581995564872319755880230778732640133618475827311280556830945027425842095283255822720234025967248535025019903452441041548430867528611936766586125109016585819976062429062342194549825839469848628615282102909112966176776453321243697465629380813815060910916610613405803768297961691307840088108557934377597967/443258066469410784554821478473651516675591396816359905867286341826785685389147294466518996332429919591238411691310164330195342302710532918831725246304362735391680355184993868960797589018949135827755953650584101506440369430035070043475502, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8+2465367414841404005816575929891574805999146616114899605374010680810604144204690768190429208855224706991025497572469763412250128983213234252887250927130677831725427027876163314121255622577744880499674013334256396478567782899843569792736628976537284185653974797587012401706396891528309108747347528418223426330841220812856004364202441418201032633451546689193069604826108382967115600595110518287525300089783777481006571972310976903856240166770398470045064654202056594110385318018518742036586681022147140900905968777522383351/783040251864423588732265591392161323435757401852876055444680926678939494359556642086819233929035877391320738064150954443907463737281485166363288225480824073861723824874400288340595260024613490089686480551020944528304983435836739968096711518991085685624842220744615721620102, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9-320998084204622756388673007698327718240335313013357797200859177162656930125817804032027891160141819624111242808277817386138963959048010340877149118319064618374543011543828478100580243364426809396904305009639689560187100894450977009722177344684651072334055906509758867840886404812625678731886472463788959615917862821582411676987186763004432931767122605495705107691651801679951034170521130844886974587964381079845207429849683459948106501499773250892083783028514221593896296036913902219539505952914271857707640360315336108322167561273790399968692272947808811416400746887357397367083041894734942817/1383284552323456743505660988225737701395192297468570418720266288662835940354876147074005754092939075619033334352636606644769470382166017076012235344659454323225370343064716046620884665855536557623665039152753826664814121717946308322547411767961750347755212294582916427342847091292636715831448546712563009344702, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -119371844055863671464916250841994507511009759491470569108577493644755390272180452429019451521584735222246375370713743316634367572764860654250030923828159172789125248054351943442726159000179642852748061562159841942233965621192369309645793888474412144340501458291170127808470618634172545120539154340970658371412923083/71018545374552022913160141760677593498913674859327751874356241069960660859738065391622108058485659890832315479417755097212213579682211575929533786755687338989963151-x2_6, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [gama, x2_6] 161
# 275, [1 = 6, -1 = 3, -5*sgm*x2_1*x4_4 = 1, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, -10*sgm*x2_3*x4_2 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, -7*sgm*x2_1*x4_6 = 1, b*x1_5*x4_0 = 1, -3*sgm*x2_2*x4_1 = 1, -sgm*x2_0*x4_1 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -7*sgm*x4_1 = 1, -sgm*x2_0*x4_4 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, x1_1*c = 1, b*c*x1_7 = 1, x1_1*x4_0 = 1, beta = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, -6*sgm*x2_1*x4_5 = 1, b*x1_1*x4_5 = 1, -sgm*x2_0*x4_2 = 1, x1_2*x4_2 = 1, -sgm*x2_3*x4_0 = 1, -21*sgm*x2_5*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, -sgm*x2_5*x4_0 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, -15*sgm*x2_4*x4_2 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, -sgm*x2_1*x4_0 = 1, -sgm*x2_0*x4_7 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, -4*sgm*x2_1*x4_3 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, x3_1*x4_1 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, -x2_2 = 1, -x2_5 = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, -sgm*x2_0*x4_0 = 1, beta*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, -10*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -sgm*x2_0*x4_3 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, -6*sgm*x2_2*x4_2 = 1, x1_1*x4_7 = 1, -3*sgm*x2_1*x4_2 = 1, -sgm*x2_0*x4_5 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -sgm*x2_2*x4_0 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -6*sgm*x2_5*x4_1 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, x3_6*x4_1 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -20*sgm*x2_3*x4_3 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, -4*sgm*x2_3*x4_1 = 1, z_aux*x3_0*c = 1, x2_1 = 1, -35*sgm*x2_4*x4_3 = 1, delta*sgm*x3_5*x4_1 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, -x2_1 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, -sgm*x2_0*x4_6 = 1, x1_1*x4_6 = 1, b*x1_2*x4_2 = 1, -x2_4 = 1, x1_1*x4_4 = 1, -x2_0 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, -sgm*x2_4*x4_0 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_3*c = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, -x4_0*sgm = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, b*x1_3*x4_4 = 1, -21*sgm*x2_2*x4_5 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -35*sgm*x2_3*x4_4 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, x1_6*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, -5*sgm*x2_4*x4_1 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -x2_3 = 1, -sgm*x2_7*x4_0 = 1, delta*sgm*x3_0*x4_3 = 1, -2*sgm*x2_1*x4_1 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, x1_4*c = 1, x2_5 = 1, delta*sgm*x3_4*x4_3 = 1, -x1_1 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.592097164
# 2
# [x2_6 = 2, gama = 43]
# [x2_6 = [1, 2], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]