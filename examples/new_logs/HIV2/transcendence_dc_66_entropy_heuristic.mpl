infolevel[Groebner]:=10;
Et_hat := [278053148019442652616-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 199544822345576027103-z_0, -c*q*w_0*y_0+h*z_0+z_1, 179611254510046812660-k_0, k_1, -w_1+1615152356527916315206945236820562205381526634675189562628379830415594401408534056, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm+x_1, -beta*v_0*x_0+a*y_0+y_1, 755383977127654839883636769068874013117736352359392881735788520575986606637183368-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+9382080920066651274464461781226633328278048002113996661360177301070616001040453959148709209153800872035349910139918823689500847249541056250280, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k_0*y_0+u*v_0+v_1, 4387866921959285211476273401129379819028755878572078725861398414429473288437170206711685842619762828315278763124556854257178887814003390543904-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+54498538193574621658397856566667706474228103391607264496019463934750657850961977522245253358427917753746380276311533183334893140735296568282821373795949223512895986916140409556480415176658151455423752712, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, -k_0*y_1-k_1*y_0+u*v_1+v_2, 25488197668734453635401566075941766301067580939793969891654509179789216000305088149481795128200236628517048451414884743392193053683137036003797188323522718703160130745477008610991230303865651552829951584-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+316570565798894426330232282144911526157833097611950292298295983915905525115010505251798059188289918561691487579500010754710019093334138999231344115529047109074475314561415023875659692000527159340933000864591463225748027267285305881908795489855166001196265548291848, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, -k_0*y_2-2*k_1*y_1-k_2*y_0+u*v_2+v_3, k_2, 148055588730206364513367662093199013645192417166270377738045264556619197412475997930371377543914444420571782820423923150973291217796678209314300633001504480129703567574474580061206954041168620687003881858870295655827573350367279026092532165970650266923482498041056-z_4, -4*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*q*c+h*z_4+z_5, -w_5+1838891949253195234853015991564776645235636954020217058533124075239352102616916878860986068193647718111110520678818113391842162920267055478456629514303174470323873749967835906135841501012388695403996027718808169537769042816955696670049024403671588614077793115238519091439269133703363745790969659073651539928269474869072348520, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, -k_0*y_3-3*k_1*y_2-3*k_2*y_1-k_3*y_0+u*v_3+v_4, k_3, 860023829034295538017503145478273329753071651205073075709003667729767315651308534480759002308630145773258653638926804684367866929802132365157801813832511554311681932461035694019868281896442688779280498074480740720265985428830142719615081981057707196927853440326237215004778970152619018829036392744167252666518441654408750816-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+10681737237619156462606416029255816906903146400472953245105189311436704780018848589168705846279918713830486113300941821363558822569561358698977053050021243722767415481838880071776987095128611631389543121598888456369902143933993132353338344491390637772506982801433309048999596965620396610985617287862669175329870759761464048972196030034850891091821738197485682759321786909198571560181224, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, -k_0*y_4-4*k_1*y_3-6*k_2*y_2-4*k_3*y_1-k_4*y_0+u*v_4+v_5, k_4, -w_7+62047968865096965740461367739291906765484835198830533561951011176360416553777109088778188740933354539829972695993857805846822565607412685568650241441950966704614561772209701947350950907823664126093116789445449499698861607999946670697424257306160924423527341698420158960333303434004484154196070293249368114687182759358217519836305077455686872329699601062268502350791953692748717730840924879425608066047089828089829579673155133063439897935169679432, 4995697851397009308222846567206048310507168674657252165527452444188307828989282915089751548590177949416985285009658481435046406343983601209571202174341526748625209457362654697172420905902799780818347310948203049036412321500174702661109002894371133717993616006107942194309817886455018201538443461475290381863479717798892194750256473260712695509849312868473842492450891680715089010297696-z_6, -k_1, -k_2, -k_3, -k_4, 4047839940045475915290206185084505077251854280701454242039263464379719789418843997272128623486912001513045885245418980200420747400330914774592-y_3, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, z_5, y_5, x_5, w_5, v_5, z_4, y_4, x_4, w_4, v_4, k_4, z_3, y_3, x_3, w_3, v_3, k_3, z_2, y_2, x_2, w_2, v_2, k_2, z_1, y_1, x_1, w_1, v_1, k_1, z_0, y_0, x_0, w_0, v_0, k_0, z_aux, w_aux, a, b, beta, c, d, h, lm, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [k, y_3] 83
# 264, [1 = 9, -1 = 2, h*z_3 = 1, a*y_4 = 1, -y_2 = 1, -30*c*w_1*x_1*y_4 = 1, -y_1 = 1, -3*c*w_0*x_2*y_1 = 1, -z_6 = 1, beta*v_0*x_2 = 1, h*z_5 = 1, -60*c*w_3*x_2*y_1 = 1, y_2 = 1, -6*c*w_2*x_0*y_2 = 1, -1238850523131833999692235920125556477868773187666833811120568482771192993504942636969435124262224497922864286069361616374697195553081269926755095129092181625603816676820701013546709562257982212444463663308468624912549712214904789937486399419798432795970920972608857389773268297997914605519006261200875892607715514323397348574636702830540004604507979310583346385652011268601801735497416 = 1, d*x_2 = 1, -beta*v_0*x_3 = 1, -10*c*w_3*x_0*y_2 = 1, -4*c*w_1*x_3*y_0 = 1, -w_5 = 1, -10*c*w_2*x_3*y_0 = 1, -5*c*q*w_1*y_4 = 1, c*q*w_5*y_0 = 1, -4*c*q*w_3*y_1 = 1, w_2 = 1, beta*v_3*x_2 = 1, -beta*x_0*v_1 = 1, c*q*w_2*y_0 = 1, h*z_0 = 1, -30*c*w_2*x_1*y_2 = 1, -beta*v_0*x_1 = 1, -c*q*w_4*y_0 = 1, -beta*v_0*x_2 = 1, -20*c*w_3*x_3*y_0 = 1, -4*q*w_1*c = 1, c*q*w_4*y_2 = 1, -10*c*w_3*x_2*y_0 = 1, -12*c*w_2*x_1*y_1 = 1, -w_6 = 1, -z_1 = 1, w_7 = 1, -10*c*w_0*x_3*y_2 = 1, x_6 = 1, -10*c*q*w_3*y_2 = 1, -w_7 = 1, w_5 = 1, w_3 = 1, -c*q*w_3*y_0 = 1, beta*v_0*x_0 = 1, -c*q*w_0 = 1, beta*v_2*x_1 = 1, a*y_1 = 1, b*w_1 = 1, q*w_1*c = 1, c*q*w_2*y_4 = 1, beta*v_1*x_3 = 1, -30*c*w_2*x_2*y_1 = 1, -3*c*w_0*x_1*y_2 = 1, -6*c*w_1*x_0*y_5 = 1, -w_1 = 1, -20*c*w_3*x_1*y_1 = 1, -15*c*w_4*x_0*y_2 = 1, -2*beta*v_1*x_1 = 1, u*v_1 = 1, x_5 = 1, -10*c*w_0*x_2 = 1, -5*c*w_4*x_0*y_1 = 1, -w_3 = 1, -c*q*w_0*y_0 = 1, -c*q*w_1*y_0 = 1, -6*beta*v_2*x_2 = 1, c*q*w_0*y_1 = 1, beta*v_0*x_3 = 1, -15*c*w_2*x_4*y_0 = 1, -60*c*w_1*x_2 = 1, -12*c*w_1*x_1*y_2 = 1, -10591982956577674406960676710223795482278585121541118595653731851017668239845049793619418663641952992570806100098981940411406122856121024134379759707789044551475609626862887049808254061461363797250430864 = 1, -3*c*q*w_1*y_2 = 1, -y_4 = 1, beta*v_2*x_3 = 1, w_6 = 1, -20*c*w_0*x_3 = 1, z_2 = 1, y_6 = 1, -6*c*w_1*x_1*y_1 = 1, -30*c*w_1*x_4*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_5 = 1, -20*c*w_1*x_1 = 1, c*q*w_1*y_4 = 1, c*q*w_1*y_1 = 1, -c*w_0*x_3*y_0 = 1, c*q*w_2*y_1 = 1, w_3*q*c = 1, a*y_0 = 1, -60*c*w_2*x_1 = 1, beta*v_2*x_2 = 1, y_4 = 1, -4*beta*v_1*x_3 = 1, -c*q*w_0*y_1 = 1, -c*w_0*x_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -10*beta*v_2*x_3 = 1, u*v_2 = 1, -c*q*w_0*y_5 = 1, u*v_3 = 1, -beta*v_0*x_4 = 1, -c*w_2*x_0*y_0 = 1, b*w_3 = 1, a*y_2 = 1, -10*c*w_2*x_0 = 1, q*w_2*c = 1, c*q*w_2*y_2 = 1, -5*beta*v_1*x_4 = 1, c*q*w_5*y_1 = 1, -12*c*w_1*x_2*y_1 = 1, -5*c*w_0*x_1*y_4 = 1, -160384182955497141382823800457813098826896606009047224230606208495201262522850319793793638262166142072894427665475312342224811258056808274298234867970937968055106695157950986057496934399896001764981724451943728043583837248557585624273089490857349933249291673230043661190995339790063375756735309317207621243238438142918919680 = 1, beta*v_1*x_2 = 1, c*q*w_1*y_5 = 1, -60*c*w_1*x_3*y_2 = 1, -z_0 = 1, -5*c*w_1*x_4*y_0 = 1, -3*c*w_1*x_0*y_2 = 1, -5*c*w_4*x_1*y_0 = 1, a = 1, -699508523128880244923001612615625014529544595235571828887627756357009846103118048 = 1, -60*c*w_3*x_1*y_2 = 1, -20*c*w_3*x_0 = 1, -c*w_0*x_2*y_0 = 1, -z_4 = 1, z_4 = 1, -c*w_3*x_0*y_0 = 1, -90*c*w_2*x_2*y_2 = 1, c*q*w_6*y_0 = 1, c*q*w_3*y_1 = 1, b*w_4 = 1, -c*w_6*x_0*y_0 = 1, beta*v_4*x_0 = 1, beta*v_2*x_0 = 1, -c*w_0*x_0*y_1 = 1, d*x_5 = 1, -c*w_5*x_0*y_0 = 1, -4*c*w_0*x_1 = 1, -w_4 = 1, -2*c*q*w_1*y_1 = 1, -5403191785151578035281168156192033510041488011599520052944548122979328446614845693097117406982629088408086692928291907877174464533530506840368 = 1, -beta*v_4*x_0 = 1, -20*c*w_1*x_3*y_1 = 1, v_1 = 1, u*v_0 = 1, d*x_1 = 1, -4*c*w_1*x_0 = 1, v_4 = 1, -c*w_4*x_0*y_0 = 1, b*w_2 = 1, w_1 = 1, c*q*w_1*y_2 = 1, -c*w_0*x_1*y_0 = 1, -c*w_0*x_0*y_0 = 1, -30*c*w_1*x_2*y_2 = 1, c*q*w_1*y_0 = 1, v_3 = 1, c*q*w_3*y_2 = 1, -10*beta*v_3*x_2 = 1, -2*c*w_1*x_0*y_1 = 1, w_4 = 1, beta*v_0*x_1 = 1, c*q*w_0*y_2 = 1, z_6 = 1, -6*c*w_0*x_2*y_2 = 1, beta*v_1*x_4 = 1, c*q*w_0*y_5 = 1, -c*q*w_0*y_2 = 1, beta*v_0*x_5 = 1, -4*beta*v_3*x_1 = 1, -beta*v_3*x_0 = 1, d*x_0 = 1, v_5 = 1, -15*c*w_0*x_2*y_4 = 1, d*x_4 = 1, -5*beta*v_4*x_1 = 1, -c*w_0*x_0*y_4 = 1, -6*c*w_5*x_1*y_0 = 1, y_1 = 1, -z_3 = 1, beta*v_4*x_1 = 1, d*x_3 = 1, -6*c*w_5*x_0*y_1 = 1, -10*q*w_2*c = 1, u*v_4 = 1, y_5 = 1, -c*q*w_2*y_0 = 1, -2*c*w_1*x_1*y_0 = 1, z_1 = 1, v_2 = 1, beta*v_1*x_1 = 1, beta*v_3*x_1 = 1, z_3 = 1, -2428543007268361829797898145935981475093566844866646368399514509232827121167119206637730541505973444561807712331375742069488412264708218901026229366441883451886798621783148723794701476063569416402563930828096268625735781494111391164762975803978446612170189733726576442663412743083375205333580904038129318008688344214247980764979670478991986220856928532633541196348240594133479116710072978645613811198700053428812385068717547776784671422075400720 = 1, -6*c*w_0*x_1*y_5 = 1, b*w_0 = 1, -w_0 = 1, c*q*w_0*y_6 = 1, x_2 = 1, a*y_5 = 1, b*w_5 = 1, -60*c*w_2*x_3*y_1 = 1, -c*w_1*x_0*y_0 = 1, x_3 = 1, -beta*v_5*x_0 = 1, b*w_6 = 1, -6*c*w_2*x_2*y_0 = 1, -c*w_0*x_4*y_0 = 1, -3*beta*v_2*x_1 = 1, x_1 = 1, -lm = 1, beta*v_0*x_4 = 1, -c*w_0*x_0 = 1, h*z_2 = 1, -z_5 = 1, c*q*w_3*y_0 = 1, -c*w_0*x_6*y_0 = 1, -5*c*w_1*x_0*y_4 = 1, -6*c*w_1*x_5*y_0 = 1, -3*c*w_1*x_2*y_0 = 1, -w_2 = 1, -15*c*w_4*x_2*y_0 = 1, -6*c*w_0*x_5*y_1 = 1, z_5 = 1, -c*q*w_0*y_4 = 1, -c*w_0*x_0*y_2 = 1, beta*v_5*x_0 = 1, -5*c*q*w_4*y_1 = 1, c*q*w_4*y_1 = 1, c*q*w_0 = 1, -z_2 = 1, -c*w_0*x_5*y_0 = 1, -c*q*w_5*y_0 = 1, -beta*v_0*x_0 = 1, -2*c*w_0*x_1*y_1 = 1, c*q*w_0*y_0 = 1, beta*x_0*v_1 = 1, -3*c*w_2*x_1*y_0 = 1, -4*c*w_0*x_3*y_1 = 1, -beta*v_2*x_0 = 1, -30*c*w_4*x_1*y_1 = 1, -81815322340113410247373836646548347083174869816821684223190817690067904156980626634559582990243776572861047930378381084997736290487009950736720452399774078602760655688831367308185001330430005011874555589058482654437908850983642811267971145276594395247095139835432 = 1, z_aux = 1, -4*c*w_3*x_0*y_1 = 1, -6*c*q*w_2*y_2 = 1, c*q*w_4*y_0 = 1, beta*v_3*x_0 = 1, c*q*w_0*y_4 = 1, -4*c*w_3*x_1*y_0 = 1, -15*c*w_0*x_4*y_2 = 1, -5*c*w_0*x_4*y_1 = 1, h*z_4 = 1, -3*beta*v_1*x_2 = 1, -y_0 = 1, h*z_1 = 1, x_4 = 1, -3*c*w_2*x_0*y_1 = 1, -15*c*w_2*x_0*y_4 = 1]
# 273, -5.531957818
# 2
# [y_3 = 20, k = 5]
# [y_3 = [4, 4, 1, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4], k = [2, 2, 2, 2, 2]]