infolevel[Groebner]:=10;
Et_hat := [7825149869257548048-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-63368340974867727071696199243738272259301410320550450671/3262529712745814377, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+2687735260518588924337371024568839409042748847272175722689918148429555546137079059406862011415965288257143134699/55749814900192883881202327666017230211627947869191665659, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-113998894708320964260839660363668361240669837441564250062033516325456206263831906027521340486340669307887689641709442385116484625296139406956728684902700677908279343591/952647833141103985499566078990323988717544484055155037790395737566220338369880933058565339953, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+36514603819611988838919862098526171346607379518614499121630763872835713198360595769928185757597344385731196284943729013014547787870752726318847814059950584011450864563670229133097431616951212205360984482936541094102580334019/16278760666975789196906380177057775826290546006288934358593144619690259495826343185470221812168773793956687059724523487932297370451, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-13456933992734523212680213606939396942093206378135868083021252030932664772719280003150484156203173275261508308822756910945972835143524373776782606857978261630682426925605949782674137998039801899225233243188232390146861088703257511888768797776795895381249462728596010503348806380561659672897828469951/278170001162882134368183867542507262269184225771511845166185549613323055159228250303443922141157287475441734081847122470049272447963881456849221005205206352428028911417, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-20657560239250642809116468307089741220546736615072898639020098979683983784274675420338819319729151527457891611474377407018703071986874843633025182814278178980236727853502415124865920900877887200871120424555791940779326837176192693019003358466133463536764649445577649315315952944221978074530805337899219499652220926360758502165349492014231385121140950017742439248944765346661/4753344012479603725266357327018968539546689650530742440597664298794175391387542160667830241390241974094018896924862826185443337345133253686176891234072604193291770008799174030639925197448327846839247615339, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-9957727203573490510009879260247376574078965648672546813048195964893915854532928190829080065707107831704230524877264004166772835148225244592311125037400896106453071022554491358144624706842600919493345113993980269160032935666734645346398890448942900184874467002699613936606330172298607590072988387388191251600291033721799439155498058371319813355054010265598561234908011377936968247309353132513219402607718697329870177403384786809800429038058187043991/81224715844702621456477902324902956319656849929772581426644191441584692825530123276316793271882307617404069775321230971571215443042703069621363370587410959894972173709870808376840118293049911932164483555356873517553670277930054783692777584513, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+5819940887626846976527167582986345550856917550448342583560424918390782840887352483034099859753969371425373101573939683524232260626273989252304912042241184797637679865445054775582292064759592639440873324218161155089993099318647319928279328800110500433500721163845871928055982672702907439240670089578417986406383197611311538562565294715541644566171416474915529669534763080214057031292397198094881612083876328726205437044301988588013855927285616104012708776569404020392819937236938349378481520562951097466281058961001937752179/1387960653958873193474917940682127291389906566650403286354200854444378001073426758276731451691368351391299868077063722611655786952811227107232016415170758852697703964629879770681156726502845139822044218566044349002104749579414833682705630468624986477595416339461786631080719227971, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+16279133765200089192707682528555346923038235815401925662141944653881890234510987171352390344660483723249479468691921612360643142609470828571604855713555277465342558603352706173874643957502885816265989229017716045331957324228996220112564000216379781813371399838016266762784697348730903647676212123412277822132153192289603278085878791906849528168286568762507332852848935676069082344427698126585338498668663987312596644540105768438357195524592762870110515242654126791606480199849140842718460274394016477318681269217191498772953030534856471151133831835734394098777193399261592251789392622463337831996049/23717347077226896854141361559617462318849781001718821331615908916584790933977585742063708321153826801353146906068143218353653363469214809213378341543338944163822122160900135800696048583025069358288101870555811083660820410446657342455919210049154159071537681193215264954692693487465544486668565195238570525989158389257, -176244840173603276900156345629273726572849624677841015139136217038929378285280981059744567694818036231769576620514309521331410042537613979403748732725822566328529953853335304229918838815245702281285959537515241412620310957033347344764148242335851956655376282068502450178206596584078591065501100045707763159175507119827935092849425949406873912037455284749101308694854318653564291515676043433301/4753344012479603725266357327018968539546689650530742440597664298794175391387542160667830241390241974094018896924862826185443337345133253686176891234072604193291770008799174030639925197448327846839247615339-x2_7, 5237625937127146571-x3_0, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_7, x3_0] 61
# 276, [1 = 6, -1 = 2, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, x4_4 = 1, -gama*sgm*x2_6*x4_0 = 1, c*z_aux = 1, x1_7*x4_2 = 1, x4_5*delta*sgm = 1, delta*sgm*x3_5*x4_0 = 1, x3_4*x4_1 = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_4*x4_0 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_1*x4_0 = 1, x2_4 = 1, delta*sgm*x3_4*x4_0 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x4_4*delta*sgm = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, -4*gama*sgm*x2_3*x4_1 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, x3_1*x4_1 = 1, -gama*x2_5 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, delta*sgm*x3_4*x4_2 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, x1_2*x4_4 = 1, delta = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_6 = 1, beta*x2_1 = 1, x4_2*delta*sgm = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*x2_1 = 1, delta*sgm*x3_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, x3_4*x4_3 = 1, b*x1_3*x4_1 = 1, x3_4*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, -alpha*x1_4 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x4_0*z_aux = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x4_8 = 1, x4_6*delta*sgm = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x3_1 = 1, delta*x4_0*sgm = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x3_6*x4_1 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x1_8*c = 1, -x1_2 = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, x2_1 = 1, delta*sgm*x3_5*x4_1 = 1, -gama*x2_3 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, x4_3 = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, -gama*x2_4 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, x3_4*x4_2 = 1, b*x1_0*x4_1 = 1, -35*gama*sgm*x2_4*x4_3 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_4 = 1, x4_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, x4_7 = 1, delta*sgm*x3_6*x4_0 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, x4_6 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x4_7*delta*sgm = 1, x1_5*x4_3 = 1, delta*x3_4 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama*sgm*x2_0*x4_7 = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x4_5 = 1, x1_4*c = 1, x2_5 = 1, x4_3*delta*sgm = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_2 = 1, -x1_1 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, x4_2 = 1, delta*x4_1*sgm = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.598868592
# 2
# [x3_0 = 19, x2_7 = 2]
# [x3_0 = [4, 2, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 4, 2, 3, 3], x2_7 = [4, 1]]