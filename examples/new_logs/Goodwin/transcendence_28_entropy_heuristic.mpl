infolevel[Groebner]:=10;
Et_hat := [1899014890301421710-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 1953190163399866286-gama_0, gama_1, -x1_1-17129664194209683532605133551966501708170589052059308849/10773764297450620635, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+973252928314736561537571094799844711655810194206098636352749710522697199424196342825647912912177886604954373/731124488939259800814146901402033209658732137563330787, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-4147284726661764399454134928256247655649962159083576188650798361117657539158824858831653986438794649411407252794768497294425708906434125219986247511118789931984763/3721143814515103534098257180568558664417022176912651197622691017424855733537055885704103705, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+61919527431296796406717649078576915237038803110567908859853345609070074753280806283834943289405019727545621992598731007753842115562826941064500663265060173384546204550801229627205323893427777627141381557298136902378365923/56817593300944892436822461613578999389129810351652458143389012249915196999302380798827582380992082976694052521826281944749142225, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-1235954067735295765752770396298917397467203821293574447307953339107178235672904209749722617803433252620463784576662211370940478975975790385953444572681969947240496895534622344386225612664544635421686892007580151471104606933495844063934529082363018232559199225691422971423124269049728007031489977/2313438426051798779835706559341736981009880749766594459886506157605002577428555557249769112964950995685596941306217728902199448907793161109137101521762820813268927, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6-905495665517528584607531666046104714450976579648088426021681472087442876664635425258214742648273505586904148747106297827409787474311426436683228673480858616466453783802660616765332697849127154653885785949159194642073665035426300914335226750949175225117879158798230231420902528344495173250382567637152418348019868332854363146195010700537447533285386789409909043802120211/19624193552988288787161033635028581532541708028312902119450397296799936362215251638179923772203765573502961210795771083563118829098768077335772339822109886126058328464370946173348473944935119541489675, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+780132379153354731835951000918144803021333257703576467698924323075337640793369872436039000587122696472166704120449462010602456358511123679792816529490949033040356930372727573171959598571211251453898041934620854387197780700915405367797918612165990407363331591624061319870901054674775633662044965410029324163748056477946923716388220104857275035551639072413434411083109388522630401581984685085749956848951988123144965779083234095056598204778277157/898916707118459372562102069537879451680066662545555207126033580143556800649631787946102058619915097973504645502034670610784838934600533783034120348389479688780485544855817908892104076765842672430481065813116133335326678228241438209033625, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8+125067129115692709446262775042761875325128893089208031299891636615113371311799481660013393917213281410362472728711143305407673193582937740875574447557010739407812652343871961610702678408365864496968853594824719465364569537604340452979784974031290761692744078626798165952517932737499658393421496754272792828159467040996609002785930098065518713080161364736696052460122958746495988515463404291628556199705194923747838898086907653129863375472820523303415886710956657983379838651182534334139572145851614595914599030309391001/183005679337108625696351544965310728956745697794553071847240249981999327727227727949781606198746541805771730505873768099403983758979905508964066142269414802316147047890245125560331957134765433226183850020883115890705884163006699555794321155080897815750427476161437522446675, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9+60130244011846919204906349577109558666949560532087357515798518568725275934152228568988716662144421277675315561070266174141843292029902310835286323785490433476718870964610784498779679908043167191601929790706370338135158036999140968921528535528501049101152327045945175904034766472660461323236409023137846755360366149666384071945160295093372427844467634775026418688058960998596670283422656275741314140963685027100637437757512602386934478997921216279905999178221558793108046066517106145279991499186319374756799600753465228277793668066711226880884673653931305170303574407697721889096371310291374977/931428863331360016896953642015758683119918939088003676608597425750483568097389216731273841320949683376474851226581382669414330891301315512874024352571211283285258732940790555472997402992232501385914949867263968570803843850935384309390284801082485731824648117545037798215287918914248923959648965992647475147625, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, 384691901343520559479843655730217157102298667943774293404-x3_2, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [gama, x3_2] 197
# 276, [1 = 6, -1 = 2, -5*sgm*x2_1*x4_4 = 1, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x4_2*delta*sgm = 1, x1_6*x4_0 = 1, -x1_8 = 1, x4_4 = 1, -10*sgm*x2_3*x4_2 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, -7*sgm*x2_1*x4_6 = 1, b*x1_5*x4_0 = 1, -3*sgm*x2_2*x4_1 = 1, -sgm*x2_0*x4_1 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -sgm*x2_0*x4_4 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, x1_1*c = 1, b*c*x1_7 = 1, x1_1*x4_0 = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, -6*sgm*x2_1*x4_5 = 1, b*x1_1*x4_5 = 1, -sgm*x2_0*x4_2 = 1, x1_2*x4_2 = 1, -sgm*x2_3*x4_0 = 1, -21*sgm*x2_5*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, -sgm*x2_5*x4_0 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x4_1*delta*sgm = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, b*x1_2*x4_3 = 1, -15*sgm*x2_4*x4_2 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, -sgm*x2_1*x4_0 = 1, -sgm*x2_0*x4_7 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, -4*sgm*x2_1*x4_3 = 1, delta*sgm*x3_0*x4_5 = 1, x3_1*x4_1 = 1, b*x1_1*x4_2 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, -x2_2 = 1, -x2_5 = 1, delta*sgm*x3_0*x4_6 = 1, delta = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, -sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, beta*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, -10*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -sgm*x2_0*x4_3 = 1, b*x1_1*x4_1 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_0*x4_5 = 1, -6*sgm*x2_2*x4_2 = 1, x1_1*x4_7 = 1, -3*sgm*x2_1*x4_2 = 1, -sgm*x2_0*x4_5 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, x4_5*delta*sgm = 1, -7*sgm*x2_6*x4_1 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -sgm*x2_2*x4_0 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, x4_3*delta*sgm = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x3_1 = 1, delta*x4_0*sgm = 1, -sgm*x2_6*x4_0 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -6*sgm*x2_5*x4_1 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, x3_6*x4_1 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -20*sgm*x2_3*x4_3 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, -4*sgm*x2_3*x4_1 = 1, z_aux*x3_0*c = 1, -sgm*x2_7*x4_0 = 1, x2_1 = 1, -35*sgm*x2_4*x4_3 = 1, delta*sgm*x3_5*x4_1 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, -x2_1 = 1, delta*sgm*x3_1*x4_5 = 1, delta*x3_3 = 1, x1_2*c = 1, b*x1_3*x4_0 = 1, -sgm*x2_0*x4_6 = 1, x1_1*x4_6 = 1, b*x1_2*x4_2 = 1, -x2_4 = 1, x1_1*x4_4 = 1, -x2_0 = 1, b*x1_1*x4_6 = 1, x4_3 = 1, x1_5*c = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, -sgm*x2_4*x4_0 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, x4_1 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_3*c = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, x4_6 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, b*x1_3*x4_4 = 1, -21*sgm*x2_2*x4_5 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -35*sgm*x2_3*x4_4 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, x1_6*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, -5*sgm*x2_4*x4_1 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, x4_4*delta*sgm = 1, -x2_3 = 1, delta*sgm*x3_0*x4_3 = 1, -2*sgm*x2_1*x4_1 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -x2_6 = 1, x4_5 = 1, x1_4*c = 1, x2_5 = 1, delta*sgm*x3_4*x4_3 = 1, -x1_1 = 1, x4_2 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.598868592
# 2
# [x3_2 = 14, gama = 43]
# [x3_2 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]