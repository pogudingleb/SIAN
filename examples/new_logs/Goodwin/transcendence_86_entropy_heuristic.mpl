infolevel[Groebner]:=10;
Et_hat := [8623263028691810230-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-367110634293176942851320010470385701774807760693317488119/6172578539325590086, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+14671282238040713482208139627776675767351293966817632237604428040767150103908591254881261214531642326626826286489/35766712519052350603501811867060915809306450418271671179, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-1172652069437746944958947654689375271643453244081462660595946613634649581556943854586787989282957011621339935048356965922572611990733644921163634230247730887445661276133/414497026249359508431734529016955052186066751649386949530629873867205291149129152310559339187, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+93728168052093330695354170673400087939048471276508860329420528344846052172747190887156805942247147607479300352133754784544082363178135728945582027507506112596007493947063797105455514637440862526264197223896660742496109108201/4803566575430800099935430214975665870723922280650725377249838448918077746605563415249602725162278666773435957217690665190630522011, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+87120182870592601523425896176294153873964595112807106135388712466317229538535231168901867243322484729553701040978088393325860151199327643611015809758955965263211203161145458042672266631872384997271878533976627517997141195874007940387533842535219224615996061571605303647872862678048774604088/55668075723936848098005669370388266492709836416674447059302418302695264762810019377913185050518689654296967956389741254841167144449396831270903803210734089108534786883, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+7038197088417512448913838112636419967126446041720961334850045366776930790485132781291116358238461250523630588261016286041336627742763627972284924054453812114025145254004262155192305066635021339832452072541249270540065872892371560905659042220211899597265065249602298062021264503919028238148648842214984318477540875190124606066796927973933218167819081986055923818519/645132029741472687212034855293203515698167042254305032177260056151745459989598582020062369912742911975158795438903901967085916112457573494398675723513492515911676672866767738455726406587944364513141270699, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+118290549196287064297781320298148237083092350445504629571130730364888683013504229982665609480083133554573582791725763888601848089562176198495745996230021422088844032766734370827102297903094187899440853732323307475101490116417364772182082857306482281997036581173406530104568725759200117956072423459883748662617200692223478645249994443239128858304639490159927097872741268381107048598664875020779583479833801581588440789917869284264760995832/7476373673527062098555719055196836327782192635298898544937626962534528899222087285535365501048196834444565083087310288380955548780042554812569537005679090726019353407537476133638349897513815946289509269632776304819774810696160658623025409747, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+181889072099249125665268312933222876887278901387260575678513736055468606109636420278049636159843438063122450374270892163125973908816318353907765701372989608523774360853261386924510191988245958790630980504286115663031422351481345779086603360211268692504149615188163290719587661653724225616200316361698365877880411407702709677337344566959449729670985675513940092953427854713699093436197671904194436563947394375248115241027925705925893129890366719220855465797734323013458525443060024626255302113997088337750536771/86642982721859453823197139426495088413489927514797965574924861132273553209316400659301684125652889719812639264776118345125026009620637566451439255142471658555352471136919270882245699928189696939401388524037794261828128780692701069991016424927513852868231860605624257605429749691, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-52618239288733037921080858573254665839287679344068263162606346651273511797982432494867903893001634275227235850689004610260641487480561576292087656339257742933753295542000758765316374640363055756545403865130285463798460208885956318286909603338947419468159927196703356891085789027059751220304389630514845569743136096803615459553926373876375427567740796239018456917680339037546866656653548065919476875972500671686845408463419264552974413958055141498265074856203540666535258175168290585815210106547804116831176437655289636896591250259760473793689443769698891936206912576526670049335862002/1004097278005490915602020853316305771069575743062252104288014661191454842249515942858936782215158366151581158405084219061666743050828085114395183452200864926314367182320514846430771197098896734630226213736765165111368750290177459769040266498490943246542463344050618266398832176796610182688708100992492795979243117923, 65259183112311480469055478784215728038585241776383701343346813173539671823598542515745575110815645232364822683344896600423255549656/35766712519052350603501811867060915809306450418271671179-x2_3, 7443146573042387912068769185531850450430664999133865191634100374800410660454842036016538978313834281895119190103312525683065857007066498530318895945853136386503780315900330563216237931432981632550316685846588686443089305072607533568848433008388842751046172378478/4803566575430800099935430214975665870723922280650725377249838448918077746605563415249602725162278666773435957217690665190630522011-x3_6, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_3, x3_6] 40
# 275, [1 = 7, -1 = 2, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, -gama*sgm*x2_7*x4_0 = 1, x1_1*x4_5 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, -gama*sgm*x2_6*x4_0 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, x4_1*delta*sgm = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_4*x4_0 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, x1_1*x4_0 = 1, x2_4 = 1, beta = 1, delta*sgm*x3_4*x4_0 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, -7*gama*sgm*x2_6*x4_1 = 1, x3_1*x4_1 = 1, -gama*x2_5 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, delta = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_6 = 1, beta*x2_1 = 1, -gama*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, x3_4*x4_4 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, -35*gama*x4_4*sgm = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, -20*gama*x4_3*sgm = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -4*gama*x4_1*sgm = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, z_aux*x3_0*c = 1, x2_1 = 1, delta*sgm*x3_5*x4_1 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, -gama*x2_4 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -gama*sgm*x2_0*x4_0 = 1, x4_1 = 1, beta*x2_4 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, x3_0*x4_7 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama*sgm*x2_0*x4_7 = 1, -gama = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*sgm*x3_0*x4_3 = 1, -10*gama*x4_2*sgm = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x1_4*c = 1, x2_5 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_2 = 1, -x1_1 = 1, x4_2 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.588688399
# 2
# [x2_3 = 8, x3_6 = 6]
# [x2_3 = [4, 1, 4, 2, 2, 4, 4, 4], x3_6 = [4, 2, 1, 4, 2, 2]]