infolevel[Groebner]:=10;
Et_hat := [4800991692883158588-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-38387790430899610205800065330873406226764254018083793719/3562541533493951327, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+479129533259992902522755786101175453408142110394630367853939756210109909389240773754793681131849336139656917640/19811508169616693782042345813280078987386682069146779641, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-5980159500327845060582619194721592871096093164550882404264300827736518018903153217716695188720969338596771159762961889688075113231679154607561697045702309767480907150/110172990901203594559631130890833148829034067626014514576772864642349313228445966755890478703, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+111234639743753757210136262980444984366612301582033364263261812355917574342533967324157798781986277368142310133161942954550005175566141941253351902694953620503891152296702703526451548289195114654074482864116489114799112750/612678642140525819851039882896839995383348872849932472163622650113955302623956902013837221494326466977465721996053571705056767849, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-1037021311987409101719368948169711410394575620329571235401092509414043452757292672154863091224223425932246906431867934607563887493266224863315366236571688330811881756348925315436268415923699759653157498493778598995331038928228020087169586959368037486385364412119783264603348176004372635739499850/3407142852931822170969275102615517532956683553286824924398595012884145625250358214385778657803074879749689908429057710783398846277296352976025877792228822764427427967, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-177939177338426003105953391599253421983211796217241730522297356450330584423934221759837928318941754071977522642722128108277080270310648480425549849883806457180578836457828649733350446114009042844504201941861847781718876511886337766632816784034388939974658501834586553442898193183515705479701173023804626940159040782523626617955396704384031596493558177864495952715855450/18947326741678401413598775634491710845821050454482690840515948256278678408940481548480560371094358165588763886648269641149161536623610161362950746109095089076601050003693827917747240204263067868229580761, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+11254240450076237338577892138822324201793400474495043143665400702650392762950066863324481181905331333352640458055961991107211556294500918246258537266984136654375737471356459204682641652778071309856574684485101464998507015754988142588552140620598554416837656710431557526781957758902041405391639849409616601476482734258523737762248592646189293608940228282922410239286179441765785114991061328679454774691999115563927178209561414334822218980413350/105367225899261458493441809871010932845228196966796179466844877397980984453090677562652080319443987400884469345691647667399492083466072511302864140502919120638584611405182850615812457725799157191792334379796430576010097065735254543754891663, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+1769535715292200151068678894220070346577896277838918171086966584333444422646456127669388026781077390718138845852204037500243121845637555577025377512510585085855709503075330645467464754238146627889034607761910275496271730297390790144556579338775774442775199486625635835286445674037452184942527259671727907549820050433916635679357353980770863279178886005555538659275704848743484795131752766064311514537847421861368940451701004261726926741925782869190370730915615534101594221792691071656641215125456495994832780647832550/585953493338159929729620040018894017704573712538976687544526060814604012949224794789247594174078173794217341076462238797887014385698336783181456262501353067876356828882114572324898721831610937670698015018774039697238318890094301496508437555569271707213703032523896938852403529, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-207437855927130034810836844264622148571850563690554541363790117768121430648911314100826048766095332409529949212766097226876229475930968096827863747561466565570309766649376360109441069797278087444303068270149073540682973340811155207855271678437798059700890109612676521414143996230827195920388397835054497593482705291368004147707147036857993344091185215309732549363052799810559538351556376412792759165905204650744289410335665414755953417305886452955860382637379027106105142354518946007959093648802307627818675331996496923737094508399117587654661938928312509411959548736768262619010152104145850/3258522689811078981877704293112488495421392075836470554502084075571244698969817220412172793893011296828866981470742581401709837908745540767183999322557185099736526172026532064986391699379989020618307162693298490333483985128973542681040934175689062800279300145568013447015047778039526059969193602711351363322069407, -5455262853528536018689834689733613574340962801405014977456269314981249217931935375432691659238882440305571040419174880753244490996809690346706028627518886933337357770479204261683962937961174929794384707852788268255020437121055145523193364079583010658119253004781075007912462766014286000542543900955717260541710400/3407142852931822170969275102615517532956683553286824924398595012884145625250358214385778657803074879749689908429057710783398846277296352976025877792228822764427427967-x2_6, -11501595592365672077321310018289016618759797194088936539693630814050596908915959273476572837343540211272615439869171281332164876230577497836759889895304107627511209051051459466906826030193526471312716686/576821941891118296123723198381325386539445380240913688883627563572509493342649040606756433-x3_5, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_6, x3_5] 34
# 275, [1 = 7, -1 = 2, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, -gama*sgm*x2_7*x4_0 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, x1_7*x4_2 = 1, x3_4*x4_1 = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_4*x4_0 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_1*x4_0 = 1, beta = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, x3_2 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, -7*gama*x4_1*sgm = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, -4*gama*sgm*x2_3*x4_1 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, x3_1*x4_1 = 1, -gama*x2_5 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, delta = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_1 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x3_6*x4_1 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, z_aux*x3_0*c = 1, x2_1 = 1, -gama*x2_3 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, x4_3 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_6*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, -gama*x2_4 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_4 = 1, x4_1 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x4_2*delta*sgm = 1, x3_0*x4_7 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x4_1*delta*sgm = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama = 1, -gama*sgm*x2_0*x4_7 = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x1_4*c = 1, x2_5 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_2 = 1, -x1_1 = 1, x4_2 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.588688399
# 2
# [x2_6 = 5, x3_5 = 8]
# [x2_6 = [4, 1, 4, 2, 2], x3_5 = [4, 2, 1, 4, 2, 2, 4, 2]]