infolevel[Groebner]:=10;
Et_hat := [2495782514384004039-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-51398989950593843297308231096498381104693248978518455359/3305260348607970240, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+31538221702175852969066289346185792892500694974983360886434094391464688471683810394772540647725268601487339773/325496345876743363533710392511979490324136309518049280, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-396710485874910334667612704762667074509041849472313389905949356924696927271872457483315843635331654132577704538602333193179519591347831954696986451585550818431209679/657113549341590842701099334168771471818197910239066299778809186884039208005349771016929280, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+6797854848546270996022254275109436870487084178784185070966476718456984677768088712979169230354405340498113195589851721998053814697303457138198274708018443969422230095611694731826944802947170749271635218299209054106422147/1473982127425969849399990227902787288295562706705185228075548861068663258500804260439052664696933674650962565276510725550899200, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+46475029584343505451758955744125754792959248890192474555857176972397988885706406970601916412825332205405769700316766669131750724886641101734833642503084660959709591182063848489690711944380241316448546230700805889657167091021077791166553343874901811874683984084676912894873415211800914538493581/4463523045750507308979146010155081579842062874518345218062839811172902074685321082623562385455029382277286626648074783261329317309068386353012703665450059091148800, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-16111049743519750010742962498055322444226121926831108658054784335690422615145157456422683875876238736203956543165165994227999713198976378457000503196419772161354194721354526203290958906529424650549285435530449046386042989828355369092607458020538048221724614678599575441328678206713535484202840209341291426435812685284138913573014254092253268903576240491813317737083931/18021962509353435664556787194158844699826205711847246866327534057676309665222212537040312165346253365912505857931315944274331152136567282403224242904375950015073254301368361438318210573716661311897600, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-490667267433649807834710228766309698617002937332570030258524896876084442406210427092078951062980816523935381009025128107831247177695564034075481513285034240596717792700413993406789680610105852634783244255417833714740647785941754313375332591519174721627217509918437208300880913457877539732544008262820164643898329182941573670755334347337794999915610629938201065868395956820583341223315358519021921642220675732756363069972444811919736559433377/20212679122614556652613609115870397686157176866836909816618093294995895736605513358875597791837779206730866193422587829499810285486940482100953549004280458758985152073581408397348568562650205011379784631890247411109776671497447800832000, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+10520250144816230376443509413724333744501558989359723363458735494459664029923782949909626302866700048350686110661541510162958862923367153524712243873362361503143712836584125975112242307564496554197771909281859861663941671590144548000947276490936948486851568147636615856951110785244730063873540711295828720387169887730697479737619342659168405459768325877419027771705935412688522134845399742766829013732103122957873473281412435501972275967536766726600581171298966954132143241674325298461016557734325221476314401513203/5440726852273570023988037516242465190195319482348165226917494679103256751234288211056427700997966370862255706308612546676807962260231592294951249866591179365473741325542015965373797972223469255948838777712519024591085963629830462999906344423025933672641685344440916377600, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+11893774979574366986579542270297820211161679143951220255886670380503894959877822397659107413045018861483292108885054190415506663856962477950409950699034492579119621825365957397512332500671724273318565760340400518132604581587595096088193019014542326038142482476953964560549432423041534249942631841217128628172908375266667058475924293757362186373878367839059093303056596086685190530169561755103534034159878944960992261917164736316110878997435263495470940651025725068485701029333134730554605576382497473852257324884325389273048940330994339420879958628484511782696736049833958377019214730314723/123567354916807108815064499820050031266190004512125358614712129664941735234906509912146306959084798707287901070042912924250832602953669368813052045115403525449169003125393830628812776979872534005314652725503748076837197085919937321527302738014189846799553486795630300882948283005206283431130154994916917248000, 108801978331152605022558966085523882039328369211126073182250790796015516066865416794546647447299465998896438692323728486279977237695700848293927063951424562720824818676333972932542936451128876115494105870282766986817477115957107995181752069/2210973191138954774099985341854180932443344060057777842113323291602994887751206390658578997045400511976443847914766088326348800-x2_5, 1342580124422626912658911621535967460070369594586729798590834273947742549055534454063556834731398221524802839852904712969689989366861478797876331776203757986300996984995587536175340228545528856303703020547204359206727783623250763601971475316037735370934272656142608029559979390311436896177038945586168697412916876166603492794825707/2231761522875253654489573005077540789921031437259172609031419905586451037342660541311781192727514691138643313324037391630664658654534193176506351832725029545574400-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_5, x3_7] 24
# 275, [1 = 7, -1 = 2, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, -gama*sgm*x2_7*x4_0 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, -gama*sgm*x2_6*x4_0 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, -21*gama*x4_2*sgm = 1, b*x1_0*x4_2 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_4*x4_0 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_1*x4_0 = 1, beta = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, -4*gama*sgm*x2_3*x4_1 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, -7*gama*sgm*x2_6*x4_1 = 1, x3_1*x4_1 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_6 = 1, beta*x2_1 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, -6*gama*x4_1*sgm = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x3_6*x4_1 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, z_aux*x3_0*c = 1, x2_1 = 1, delta*sgm*x3_5*x4_1 = 1, -gama*x2_3 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, -gama*x2_4 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -gama*sgm*x2_0*x4_0 = 1, x4_1 = 1, beta*x2_4 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama*sgm*x2_0*x4_7 = 1, -gama = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x1_4*c = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_2 = 1, -x1_1 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.588688399
# 2
# [x2_5 = 6, x3_7 = 3]
# [x2_5 = [4, 1, 4, 2, 2, 4], x3_7 = [4, 2, 1]]