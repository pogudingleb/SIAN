infolevel[Groebner]:=10;
Et_hat := [2352105420011062993-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 1883955928541425154-gama_0, gama_1, -x1_1-117569947271115666509774748534764509568432728450950357071/7342211876847378537, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+4785062338440926378169348773758719083804451071508972511709773732983769238196106655795156109595902550104317078501/43894040671957861437074704255673316003251261470728583014, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-48687658101023468073805455228532297648199172117004107638133867041446401494481148095261212481801161018996257670017231876106936112922655825657701280973221917708435052350/65603078432911474369630733435295741912819072913643365811780942353738229005623487934514090177, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+660507216161057270113342432272117466099704052716680819163567319145205615965958849763540527043299384768715982517118447570456991984937603939624131045071867213326243960865792569035768884140336106654603361302698233174575977679/130731912092240998988554715579694379995948980013402091629302186212109911723858325129163044747935967646217339091467031091661327698, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5-242510238328035216537667775577003144146179196430426428734185092163256843911127860213116791668944886464132141587034819916735892737616080501776957250878709915612890240776431206828557454547575140043536106244598297688278595182035269293849537788205558520122410059601106439574835410858006581513935/781556288860952525539368549771476250130305986887130727342838817222751715356688685339158564351237699872575533789561680647871176715663255930807557419485423584849310156, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6+799320969086247242798973944532139761261156154112647605646807042524418699288621949571664523479337949167481685212703994909615099277275797214660920516422997585813314251381348048453258614276618817438704433011136529336301321042197706949509044569397679923807127519554018826920075400023382434730616755589450992786731136882538361405484176032416399981937231129170083588013/2336194058827499276213590533229844861582604689608689850588173383256878043965308352634219064680655476861977534360716518330327126849087650028595723099297323019205795306407560103396738456265658475287360516, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7+3181405339070739613961072643522670628740963798204380097979190433499179052260539969904883551553709870021770733672649275160694326010052572478103387151151482551665610812487739602276310218329371387054993979277748572835633464450432226004066015461924650877763484767748160395740602176862511535339327554939130202957643401582441620964563292780975394967329079880842412696147621694299197502673023114439848123759411675598212393829598079800407908737/1163874941636248624313513380026281211954191528687547572542827138574589826313433879848421072205342304805491196013218673041070125899093729583918906719000460116608482689774486114275210556212354697352526680062032864251681300366170374178273746, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8+2727136068783031473790133126752890207414906508486876334199109323137194372810483364243600665005680940409226440053311490391047600755425502997367567824983614456408975614926059236596613969142842167960813169102309984432716645078660705220628279419775379216413755217225353564055593319722272481519473019156252779779758823868393987199753831804530015996619737018538095535300098877299916581494168620769630171198092982341475534434781655050182570341512614112304834930001836904050242216254421594687232570927190522155020063/13916017375186412314732299431342076941769745968185777079965420515356805720862318043860258876490120622223654093473283817638273739192990304276988358360512624026875035068806474063216806077247708394402495423042738840784518620884922822173137552921769625682013578114848205846768024, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9-791144448903141844962950556824192901582490419479108989136377555008312560052447336974737483448196329046069781568078385472937621939543590443725634888575192027942368322577004831734631445786298760155373201468819828476682892954874826565096372774097418638755957867779973345171234923730953018404360955203502448344690837752803882519843950347324787917729984004383485122156850285962293091718864158522125728151235008278969032844523286620363937460249074242803585405143669670386237974881404673242922741947309952297864115359314082342694932624310065740419495563475437415767289885286833675247505387/20798576876549659861976471957390586772100596261188955896794210665730580519459449690911289016244049153659283020477811099239830384328834177959834967615997769082775572763998885184463613194589650977805027227240164979458569194731916274163032561438879554246768733661199574872403826400608282149174798018412402027881732, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -354594903438412218037596829322905313612524456156754691391113863665181322813841901169652013946940125105174943816646741259117079370512637974604507380260544079123697711296963832124555657764527806808971378402382337495988184285917384107081521703856883848429016690590627717771494823965435563030562879083634464993143892681179983730027/65129690738412710461614045814289687510858832240594227278569901435229309613057390444929880362603141656047961149130140053989264726305271327567296451623785298737442513-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [gama, x3_7] 165
# 276, [1 = 6, -1 = 2, -5*sgm*x2_1*x4_4 = 1, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, -10*sgm*x2_3*x4_2 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, -7*sgm*x2_1*x4_6 = 1, b*x1_5*x4_0 = 1, -3*sgm*x2_2*x4_1 = 1, -sgm*x2_0*x4_1 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -sgm*x2_0*x4_4 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, x1_1*c = 1, b*c*x1_7 = 1, x1_1*x4_0 = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, -6*sgm*x2_1*x4_5 = 1, b*x1_1*x4_5 = 1, -sgm*x2_0*x4_2 = 1, x1_2*x4_2 = 1, -sgm*x2_3*x4_0 = 1, -21*sgm*x2_5*x4_2 = 1, b*x1_8*x4_0 = 1, -sgm*x2_5*x4_0 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, -15*sgm*x2_4*x4_2 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, -sgm*x2_1*x4_0 = 1, -sgm*x2_0*x4_7 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, -4*sgm*x2_1*x4_3 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, x3_1*x4_1 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, -x2_2 = 1, -x2_5 = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, -sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, beta*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, -10*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -sgm*x2_0*x4_3 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, -6*sgm*x2_2*x4_2 = 1, x1_1*x4_7 = 1, -3*sgm*x2_1*x4_2 = 1, -sgm*x2_0*x4_5 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, -7*sgm*x2_6*x4_1 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, b*c*x1_8 = 1, -sgm*x2_2*x4_0 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, -sgm*x2_6*x4_0 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -6*sgm*x2_5*x4_1 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, x3_6*x4_1 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -20*sgm*x2_3*x4_3 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, -4*sgm*x2_3*x4_1 = 1, z_aux*x3_0*c = 1, -sgm*x2_7*x4_0 = 1, x2_1 = 1, -35*sgm*x2_4*x4_3 = 1, delta*sgm*x3_5*x4_1 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, -x2_1 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, -sgm*x2_0*x4_6 = 1, x1_1*x4_6 = 1, b*x1_2*x4_2 = 1, -x2_4 = 1, x1_1*x4_4 = 1, -x2_0 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, -sgm*x2_4*x4_0 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, x4_1 = 1, beta*x2_4 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_3*c = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, b*x1_3*x4_4 = 1, -21*sgm*x2_2*x4_5 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -35*sgm*x2_3*x4_4 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, x1_6*x4_3 = 1, -15*sgm*x2_2*x4_4 = 1, -5*sgm*x2_4*x4_1 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -x2_3 = 1, delta*sgm*x3_0*x4_3 = 1, -2*sgm*x2_1*x4_1 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -x2_6 = 1, x1_4*c = 1, x2_5 = 1, delta*sgm*x3_4*x4_3 = 1, -x1_1 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.598868592
# 2
# [x3_7 = 3, gama = 43]
# [x3_7 = [4, 2, 1], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]