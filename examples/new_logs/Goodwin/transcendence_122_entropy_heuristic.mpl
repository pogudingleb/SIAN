infolevel[Groebner]:=10;
Et_hat := [6830665250930799689-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-557243854299266675081450062761624234230400367800281345343/10731798083615139929, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+141861846814664915760866972003932702806672614986054387870730043700852106105171396777295577026172973628222920841024/359404096394909323725769269925231867088609792070285932833, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-36114859637123032929155617912738373009197604847240921353797863227297571763957875132109715786914070613253555078660704556784862230415218806521180860718275915053314067771392/12036315210090899272043203873975599942740203026126053257383593712892592208982768984371141427241, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+9261184488035865187587922884232642559254975241277065631924946686006362799562939395467351809551475768688801826217290101665582990500807704349357418791656364772463336439734318587338264381215616884582188825877222831846455119840576/403091910442447419767848403292501753989590655761821313749408866817376306832567455297189217723025914811814398215745503663295466598257, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+46214081135014510308076711293760416706360611717299699729534407390276708037606699799532077213815663658828479596414937026722981042668675343235447797101637379569348596902813247671601799040978161954985170228387181953728067533409483344975958781052281126400692441441718765432496079361634165104972184630912/13499404545996013772316245568557182892364247768776463494521437237612194645582979745701542137343533295950807447776086105162927509158717955500896832042474749677936430909689, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-56325600065399104145183132508381071132712685087755571122874827793770384889446853223840068518124189674369510964929769398108465602514538728333136507057871374210670879526368304817042708885421578148834700501569231803622110373412803024046306458344638805851739215646726291064335421023916117308299627130570414717744123385247550907693514226888480303578491697563391306917987876222016/452090251318692249411841519763917608364392219866576022734602082397593479014325992891920942825926901243273564698984681968251357046086558107103441959418524841408444557610990636478592366835383129812085462556353, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-96181428014437515444458073073749614870748353331961991428033090578503468018302416140118125565048001496098476110599634253319991640312466548651005071782497339066709901932694159269871757912238262613647809621221712430969357824863621141150398110174793116480446736497007686167581747931854979908161401055600635957479914686509388424401139781786668681379973171662243217849431676366510599761222089860393061055169900226481720115443451604864435562210956169974272/15140341534399014521306037635796816401874641033411065359601961402526659575148423053564102608895276705188165903076263068529242281372489115849435563073440507770650891124866734204086422810872677655475236823071594654304932617230711655829158792942281, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+83820130994970480413120641732485159423874316506759697466242513544750410307301824085022393299205281015377699690308686035495260713465076441748096880295167292609803639621127377835403905489827786332943386769884122599116106872085088287545601318815729651117181357041100577965657166550634501791925806395776237759335022087498694669535129622483686372151045389546645312902788350580363191922407346169670654801230375232912143532343425979508968595612884205795772295247901536479035865808754697747546971617570745207665733449011893632798016/507044646748325272758205784908236023335440453218454919220286793072825440745481987794093186436394100152149912311082724786858357272049589505043692347669996479726899876775764063908396602535733829690960915709022720346789655027862113020635943065062960617514328546522135255097104259864337, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+332906148387366669758110757975890618617973731438516817917744446362856452357592891344952820297942308335955933913668406977565646891133534881554351305893668167822570615234134368764401062654086988063353061987473491561995859155930602266798043993920874370235980777317260224186764817297975162924515009706648428925793099419812062700100134678181385010451174094190593829141794516443247631582260853831304015424634659132052323143980021783540252717945462908660443989695030022989953376133542842238861092616614287264212816766238896748854808088764788559270859694324798309598238418481021147921030920985562813575748992/16980744668937162358438676756367937094714725546202036470732047967944102423814492252734807184259151538416429707726502238505321556803965872838927537174738002709188607878375315794157562952844277033768481893068535335209326291331191183491622218213169637880364593370790864109475308417929662644783482879083433125653343802813849, -301550823133027178889236502679421833715767717035148535190909785200274895836445161903026759422385668489739997197502922989581543823903294133476268663099153467536672461961382410748601391790488208586181457375916233880197921320306722771021585586404246596878122112599598686908486965245915490798625353542008113067091453683475971978706440592824733826476564185252349047684887539718739290196465354189397/452090251318692249411841519763917608364392219866576022734602082397593479014325992891920942825926901243273564698984681968251357046086558107103441959418524841408444557610990636478592366835383129812085462556353-x2_7, 7654957674862513726588121734480476691962385429185095248034526329454852122069361337825754777139665678308021977101240233470898766407966741603476382394014/359404096394909323725769269925231867088609792070285932833-x3_4, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_7, x3_4] 32
# 275, [1 = 7, -1 = 2, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, x4_4 = 1, -gama*sgm*x2_6*x4_0 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_4*x4_0 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_1*x4_0 = 1, x2_4 = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, -4*gama*sgm*x2_3*x4_1 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, -7*gama*sgm*x2_6*x4_1 = 1, x3_1*x4_1 = 1, -gama*x2_5 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, delta*sgm*x3_0*x4_6 = 1, delta = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_6 = 1, x4_2*delta*sgm = 1, beta*x2_1 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*x2_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x3_6*x4_1 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, z_aux*x3_0*c = 1, x2_1 = 1, delta*sgm*x3_5*x4_1 = 1, -gama*x2_3 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, x4_3 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, -gama*x2_4 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, b*x1_0*x4_1 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_4 = 1, x4_1 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama*sgm*x2_0*x4_7 = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x1_4*c = 1, x4_3*delta*sgm = 1, x2_5 = 1, -gama*x2_2 = 1, -x1_1 = 1, x4_2 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, delta*x4_1*sgm = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.588688399
# 2
# [x3_4 = 10, x2_7 = 2]
# [x3_4 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2], x2_7 = [4, 1]]