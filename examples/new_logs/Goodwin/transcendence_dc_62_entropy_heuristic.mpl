infolevel[Groebner]:=10;
Et_hat := [2449016904652537578-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-57403962009024725173356518593357548742369872079232146439/4681597898080076849, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+460135229426888257540794125794391807739834329902688793039534655626240497565753510335716771558937202820811074050/7495174538985559574056549101573552509388029739201805193, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-3688324323788817574851760218034016690526028610354321895739972313672492479498250877126084655534145795628539612507124360524791618466313194138004307699341442293612298875/11999672460741629360579614261484179941496599205788352518529684937575265312937596326957069601, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0)*gama)*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+59129293554305476539352984088416008271010950784412980735601865853581059336114442291712971581206535384637168059806298562362123963110312931843820398429363316222692976249927700648729083985192536075968500797591366430865477225/38422624694359633233652487698966063492232394552340314359589529448543687364588316243194979208094836598308541274296631249455164914, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+2318937694515796574303428163332229141647045972018145785659563349293411630545874487335986858657057429953175722035134321675949548692541947906685024054442804422972915781913994339644972828552531120333276073230834161594167166122846335049636080887951829986673056507951771152480765241395574415/61514099373691295520331418120602400931717676089438126590192007184440079955149482903478643899060725297845329240101143993903618894555267621718908397743058531024612898, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-12978855861981165401470076986894218338923497596296221029714931124736904340486687697387475241637603927515992188843567014616312329790766517510605977807585910378515426362215683083792001447302574098890696871143410818166566247341225570054805219267404991484471280863837677842850751447304185706332498459739742340805638210873567753646218613215000137813844067900455/49241618081232350169067842081509132762380420510239839028362378873232992494547070636546186697147945460216010564170781025758005770124788008554479234673473433074885198213503092249109766753464373656812393, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-1030672821698198114074938469053900360587030856569351461595504867629636928435985970577309439932580203025058738767262003637515327191607355565755017080308129281483071120435138769025323481898962997746431718588495228001205930195446551024108686119524785432182879972348482503059272697509072527422851673321243914501030395744919123151716616465168448222297774806771628977587262045865409271486750622751010062130882838445254582751063629867005/315340642349708098353082233640070693852706320457538646331219116441025909768665976293071776488840888678549902482476037038217389909574535158157372526447276675399886675337501562957579567157145261431197267636780301812814213671541781991680004, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+173320130681629787080506704976015705312389883570685951050870050056383707372459606628823289110078541152824511538422954384213129715974405892009084608557117011221497746802530545589400937279307081176615444302083760186932493849758854025785209274505384364696882472942056690807125899050815247600953212161316062430146856687524321627131272725597575009553956783331355702059000828895155730986891978492260851009810510009941700496044319117044212663607517948466562350636176420112786888790135288766626448602338128295/252428039005247384847033693682983202222693297790635849961954326705587788106834474518365526337059908826952742402724135725469228720948648326988794063415046921364450941259687090819236610224772398049734882779047611523893617923661496499978208714071973362438362320216431908270514, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+76962125484968597750297267636153272742855446668587308468142053307146917936705137816247869233824528145759855951964592733500833168521678006110234887445058258626681907903222816528380165159833976277647818066823246248324875795398012439434556375850753645784339182200695702859187313551642358517488436642030088838914050271968950357161866926615716996607006593883623433015572198122233590056947831355035615900536034086883458931309964584636535739008229295479014107579847290677468808497784050331760032209171159191584011693825117955908392129710950039173220631789596441554333999378492650/202066928009140838875858519723320034237981821295270030390455950267619867872099066542599852753183863578663089653736534850300486599655711440430581529113669254416592316154832837474636380435122018612805304681125537031348263023714217009628649236961071758397486707673992476181947117469732112790788905723111332376049, -19590057372632498535279180453454374450-x2_1, 9508058540979713484143344526165405480361930811305791704123116995976554815584979296400196118578405904961817052670817761620204126442472324387154021003497720979031027322768788380539541949229377765260998828480243790584451879823826666954777274025580355519318104602472229377381621517091811730067576585091094949991343617717933430/30757049686845647760165709060301200465858838044719063295096003592220039977574741451739321949530362648922664620050571996951809447277633810859454198871529265512306449-x3_7, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_1, x3_7] 40
# 274, [1 = 8, -1 = 2, x3_1*x4_5 = 1, delta*x3_2 = 1, x2_2 = 1, x3_1 = 1, x2_5 = 1, x1_3*x4_1 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_5 = 1, -gama = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, delta*x3_0 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, delta*sgm*x3_0*x4_0 = 1, -gama*sgm*x2_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x3_0*x4_3 = 1, -gama*x2_2 = 1, -21*gama*sgm*x2_5*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x3_4 = 1, delta*sgm*x3_4*x4_1 = 1, x3_1*x4_4 = 1, -2*gama*x4_1*sgm = 1, b*c*x1_7 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_1*x4_3 = 1, -7*gama*x4_6*sgm = 1, x3_5*x4_2 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, x3_0*x4_5 = 1, b*x1_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_7*x4_0 = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -gama*x2_6 = 1, x1_3*x4_5 = 1, -gama*sgm*x2_0*x4_3 = 1, -alpha*x1_2 = 1, x1_6*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, -x1_4 = 1, -gama*x2_0 = 1, x1_1*c = 1, x1_1*x4_6 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, -4*gama*x4_3*sgm = 1, delta*x4_0*sgm = 1, -gama*x2_5 = 1, b*x1_3*x4_2 = 1, -gama*sgm*x2_0*x4_6 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, -gama*sgm*x2_7*x4_0 = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, x4_1 = 1, delta*sgm*x3_4*x4_3 = 1, x3_1*x4_7 = 1, delta*sgm*x3_0*x4_7 = 1, x1_1*x4_7 = 1, x3_0*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, b*x1_3*x4_0 = 1, beta*x2_5 = 1, delta*sgm*x3_3*x4_2 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, -4*gama*sgm*x2_3*x4_1 = 1, delta*x3_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_3*x4_2 = 1, x1_3*x4_4 = 1, delta*sgm*x3_1*x4_3 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_4*x4_0 = 1, -x1_0 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, x2_3 = 1, x1_9*x4_0 = 1, -15*gama*sgm*x2_2*x4_4 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_4*x4_0 = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, -20*gama*sgm*x2_3*x4_3 = 1, x3_2*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, beta*x2_4 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_2 = 1, x2_4 = 1, delta*x3_5 = 1, x1_5*x4_4 = 1, -gama*sgm*x2_0*x4_7 = 1, x1_5*x4_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, -gama*sgm*x2_0*x4_2 = 1, -alpha*x1_5 = 1, delta*sgm*x3_5*x4_1 = 1, x3_3*x4_5 = 1, b*c*x1_5 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -gama*sgm*x4_0 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x3_1*x4_6 = 1, b*x1_0*x4_7 = 1, b*x1_1*x4_2 = 1, x3_0*x4_6 = 1, x1_3*c = 1, x1_4*x4_2 = 1, x1_5*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x3_3 = 1, beta*x2_0 = 1, x3_0*x4_8 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, x3_0*x4_1 = 1, z_aux*x3_0*c = 1, beta = 1, -5*gama*x4_4*sgm = 1, b*x1_5*x4_1 = 1, b*x1_1*x4_7 = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, x1_3*x4_0 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_1*x4_1 = 1, x3_5*x4_1 = 1, b*x1_1*x4_4 = 1, x3_1*x4_1 = 1, x3_3*x4_1 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_0*x4_6 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_5*c = 1, x1_6*x4_2 = 1, -10*gama*sgm*x2_3*x4_2 = 1, b*x1_1*x4_5 = 1, delta*sgm*x3_1*x4_6 = 1, b*c*x1_1 = 1, x3_5 = 1, b*c*x1_6 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, delta*sgm*x3_1*x4_0 = 1, -6*gama*x4_5*sgm = 1, x1_8*x4_1 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, b*x1_4*x4_1 = 1, -3*x4_2*gama*sgm = 1, b*x1_2*x4_1 = 1, -alpha*x1_3 = 1, -gama*x2_4 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, b*x1_4*x4_0 = 1, -gama*sgm*x2_2*x4_0 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_5 = 1, x1_4*x4_0 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, delta*sgm*x3_5*x4_0 = 1, x3_0*x4_4 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, x3_5*x4_3 = 1, x1_7*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, x3_0*x4_2 = 1, -gama*x2_3 = 1, x3_6*x4_2 = 1, x3_2 = 1]
# 282, -5.577999884
# 2
# [x2_1 = 10, x3_7 = 3]
# [x2_1 = [4, 1, 4, 2, 2, 4, 4, 4, 4, 4], x3_7 = [4, 2, 1]]