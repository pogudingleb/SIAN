infolevel[Groebner]:=10;
Et_hat := [3247208199439219489-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-192497444277759529396657476619994264063750839924939501049/8863729180315075315, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+536181667300504684097567791141498820187634035461917370980015898945360699396846016118531361482168389840283299946/3691519072221325695329710303802127330981100820827126275, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-1493478427351592413178173572899753187727917053891820715620351091459110074452401970363828544362168376133493847968566443023554242288577492398783016731510833504368759908/1537424348527917543393903138924112961423582154538405060437266617191838878706876376108665875, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+1386777128996844692890130753093890547343101414475617490100574441533258126083263690311149769158858285375204435452396531470190675495145955939163613163502482140651094365682976639874381534574162044057213613297062756708895592/213432788435988445504880863311299567446405643183084291342130300866727794324050750262623677901121517217787676833294469326870625, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-185494796997871320480360813592860529015507923705237411759814145137466153469785293589195667934515486505652248175316074576985415831366461784748063851851646142379919711306573839999124607854823991066521699608246020842953535999523577324709659276849557109143657778737605683305185592535020723632/88889359446880541439532153990061165410201504855441963239956680328393888211270770992440703959831679763192746313470400941693830130459310206231791879194700523940625, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-10300312744618115894040880839382480350690657601681600953982678592564852676556558222284588221204459692384294443318269405910873065369810609894852801006600256302889320451693488392453363905241381518817873839072543195343516575503631796653975714790717045955648559727370597238549058929370711245681498261145903410027699208785222597989797371520594548534354092940204832/7404034104391087612813410580119755159527634345234274736920214241826788543650653424094932271672168947179438576208946422073292343503401915538906336655506711285127114879808915311539186968182107178125, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+118778696783381543936769213506607355047617400814389130900215330890707076896574251182047844145169747547647779260091163264904029160171383184714185626332182464274945524599542796455116886292881420345376114993754461039640176016579199328472614776069078986928151308398279297451923710770891722789348578873600578998294421141008882787830365081652920228683467261106335506933472494494221622968314428155109036034042960677526056195718838407667136/1027864327821038632535118138360174254605043142610999102805517982759453378515659808266339686582227982191536093769953469800935562040350254778857713439537469634081821940626870048837391880801557205989996888515643622507902144114269109375, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-572028556396728462679753092410003319931670007296730394350815122117811938027490719769812418262241712833652573665357851319194412818008256186227167214685371831576433680840854808215481450819009173358710212160215193392129327118228641786715727714357364773162201968919546026928955946307979370967919004002911908314757956928131220234255924715197907157422378625666156035020656928985051310435988319526799467715917535481527235976417836612230811658726041072414481223375376263565323431019515818403523642405562893662848/428079501597818427266164436671896985740826667780587983487710240448369229694292076674622536719148309112341240191820954934306989097183696334378438957553020784471486086413921457594619516690226562877639564439012208595160738970556681809219734181586914520165129945888359375, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+1524979348972454684240388409397767141049799485915230552187709357687165696636592015586343273169018244260681839558989551533068945045465691096369037613365720216442415352353582264101941433662754533841707452159869386254449880437291881963411122149145839773651286771687760993423405516909468934564180778727777993397603115135418365041323079045176462852596687782180834866728447347016243804840909360828748658891444534849002186782800592412577705913566409847432474898854888685392328503514893224917998466405635407989396528819260515969305276347652803530700750812910148598919658941692436841216/178284287846345640835665266828362257328020936745453979226008914693928094149145444279765766901775709409926464995362623668176246540391390847216710762275692553774555254039285123086292139088905687998150749416052005657544194815413571873933650891066770339629587018523028836785121470061219360711547232799609375, -2395037332831479271603365885979979923524630642861754563668056972421190830882474915286885728003502743797960079674611677868176086977239598309997726448332794699367671061576938805604893288/102494956568527836226260209261607530761572143635893670695817774479455925247125091740577725-x2_4, 2749023936727312229616540527845849380781773333999359451026064739028766867032174425533154768268446418825048255758388511072402518056717232080457090496/246101271481421713021980686920141822065406721388475085-x3_4, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279);
# [x2_4, x3_4] 48
# 275, [1 = 7, -1 = 2, x3_5*x4_3 = 1, -16953024780643622587392328659282669733570887695561802252882709385663669317199715519031492773427613034142883121153025106781108389333448928835436567072120853953844784881239863460182687590296534939053434470892731630188283515996679236064033760141356571641441377570471387029156756044200851786401349971499669980198703877723197011223711469335369557117634278409950707499158084814959121939144255104946691312798903014505840183252605710136919387709228036353411398596317037227852936180683573242310574532225856307572237228532310319385250425882547625645995974319745261693138484391284434349751603834635481063/231114970048606662748159931787522807013080339730424114498798412302971645000652402589746243646244257678711027726231794992879021344300220425274567522885843550422872869571010229858117268621009941304986343140760276638784240202637566765165609682391700649892149071714579006081022807744208829727259212279163601515631185000 = 1, -gama*sgm*x2_7*x4_0 = 1, x1_1*x4_5 = 1, x3_6 = 1, b*x1_0*x4_4 = 1, -x1_4 = 1, x1_3*x4_1 = 1, x3_5 = 1, x1_4*x4_1 = 1, x1_6*x4_0 = 1, -x1_8 = 1, x4_4 = 1, -gama*sgm*x2_6*x4_0 = 1, x1_7*x4_2 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_5*x4_0 = 1, -x1_9 = 1, b*x1_0*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, beta*x2_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_4*x4_2 = 1, b*c*x1_1 = 1, b*x1_1*x4_7 = 1, x1_3*x4_4 = 1, x1_4*x4_2 = 1, delta*sgm*x3_0*x4_0 = 1, x3_3 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x1_1*c = 1, b*c*x1_7 = 1, -gama*sgm*x2_3*x4_0 = 1, x1_1*x4_0 = 1, delta*sgm*x4_3 = 1, beta = 1, x3_0*x4_1 = 1, x1_1*x4_1 = 1, x1_8*x4_1 = 1, b*x1_0*x4_8 = 1, x1_4*x4_3 = 1, x3_3*x4_2 = 1, x3_1*x4_5 = 1, b*x1_5*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_5*x4_2 = 1, x2_6 = 1, x1_6*c = 1, b*x1_7*x4_0 = 1, b*x1_1*x4_5 = 1, x1_2*x4_2 = 1, b*x1_8*x4_0 = 1, x3_7*x4_1 = 1, delta*sgm*x3_2*x4_1 = 1, b*c*x1_2 = 1, x3_6*x4_2 = 1, x3_3*x4_1 = 1, x2_2 = 1, b*x1_4*x4_4 = 1, delta*x3_5 = 1, x3_2 = 1, -gama*sgm*x2_5*x4_0 = 1, b*x1_2*x4_3 = 1, x3_2*x4_5 = 1, delta*sgm*x3_2*x4_0 = 1, x1_4*x4_0 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_3*x4_3 = 1, x1_1*x4_2 = 1, -4*gama*sgm*x2_3*x4_1 = 1, x3_0*x4_8 = 1, x1_9*x4_0 = 1, b*x1_4*x4_0 = 1, -3696690560158270646842702889116997070449711891044102010828652700383901229930505647385455394476864342027995045173521789384491131607607460223183074178530115258646278617658022020031311292149664954312497231402781703261905255025218350581414987452427351309500900317046887668857634056194817739833525723711677395172344540038882595416750139755260480285698156333100493613474121086370620956920303334455033056209358386748273425312323408711263563946502193603859857175291957714833584936778945474544792763041403176961219020983473903837/151378130238239209601614127747050005610472817776301062288796829094847842349137378798686731149751968107566328837917468214725233118892395702511617908865682169379553564835720449449807889602285753499692788537041452202727354983528783501369621966462376751509589618549268876690383817500 = 1, delta*x3_2 = 1, delta*sgm*x3_0*x4_5 = 1, -7*gama*sgm*x2_6*x4_1 = 1, x3_1*x4_1 = 1, -gama*x2_5 = 1, b*x1_1*x4_2 = 1, x3_2*x4_4 = 1, x1_5*x4_0 = 1, b*c*x1_4 = 1, x1_3*x4_6 = 1, delta*x3_0 = 1, x1_2*x4_4 = 1, delta*sgm*x3_0*x4_6 = 1, delta = 1, -5*gama*x4_1*sgm = 1, x1_5*x4_2 = 1, b*x1_2*x4_1 = 1, beta*x2_6 = 1, beta*x2_1 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*x2_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x1_3*x4_0 = 1, x1_9*c = 1, b*x1_3*x4_1 = 1, x3_0*x4_4 = 1, delta*sgm*x3_6*x4_1 = 1, b*x1_4*x4_1 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*sgm*x2_0*x4_1 = 1, b*x1_1*x4_1 = 1, delta*sgm*x3_2*x4_5 = 1, x2_7 = 1, -alpha*x1_4 = 1, delta*sgm*x4_1 = 1, b*x1_5*x4_3 = 1, x3_2*x4_1 = 1, x3_0*x4_5 = 1, x1_1*x4_7 = 1, b*x1_0*x4_6 = 1, x1_1*x4_3 = 1, b*x1_0*x4_7 = 1, x3_3*x4_4 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_8 = 1, -alpha*x1_0 = 1, b*x1_6*x4_0 = 1, delta*x3_6 = 1, b*x1_6*x4_1 = 1, delta*sgm*x3_3*x4_0 = 1, x3_1 = 1, x3_2*x4_3 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_3*x4_1 = 1, delta*x4_0*sgm = 1, delta*x3_1 = 1, b*c*x1_5 = 1, b*x1_2*x4_0 = 1, -35*gama*x4_3*sgm = 1, delta*sgm*x3_1*x4_4 = 1, x1_3*x4_3 = 1, x1_2*x4_7 = 1, beta*x2_0 = 1, beta*x2_5 = 1, -4450714391761638746958026096852287791946222252597595940124437484399624497584321575678719708671729158299220887506356198377280161585946988503720434421763355305407467717337457157373111383965657322640149191299933830063549835516671088640534480376601888645882990871186987647648807026013693644677134019/98027983535522127929044014599770895396030734414558153746123274772181271066289530007865226247847307434414464498515781090868864894809960762585444124208018751428080363825 = 1, -x1_6 = 1, x3_3*x4_3 = 1, x1_6*x4_1 = 1, x3_1*x4_4 = 1, -gama*x2_0 = 1, x3_6*x4_1 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x1_8*c = 1, -x1_2 = 1, x3_0*x4_3 = 1, -x1_5 = 1, -gama*sgm*x2_0*x4_6 = 1, b*x1_0*x4_3 = 1, -alpha*x1_3 = 1, z_aux*x3_0*c = 1, x2_1 = 1, delta*sgm*x3_5*x4_1 = 1, -gama*x2_3 = 1, -1088856666242691992898892047617511110126853435302653987826/15839546878037681617 = 1, x1_4*x4_5 = 1, b*x1_0*x4_5 = 1, delta*sgm*x3_1*x4_5 = 1, x3_2*x4_6 = 1, delta*x3_3 = 1, x1_2*c = 1, -15*gama*x4_2*sgm = 1, delta*sgm*x3_2*x4_3 = 1, b*x1_3*x4_0 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_2*x4_2 = 1, x1_1*x4_4 = 1, x4_3 = 1, b*x1_1*x4_6 = 1, x1_5*c = 1, -6*gama*sgm*x2_1*x4_5 = 1, delta*sgm*x3_2*x4_2 = 1, delta*sgm*x4_2 = 1, b*x1_5*x4_2 = 1, x3_7 = 1, x3_3*x4_5 = 1, delta*sgm*x3_3*x4_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_6*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_6 = 1, b*c*x1_0 = 1, x1_7*x4_0 = 1, b*x1_3*x4_3 = 1, x3_1*x4_3 = 1, b*x1_2*x4_6 = 1, delta*sgm*x3_1*x4_2 = 1, x3_1*x4_2 = 1, x1_7*x4_1 = 1, b*c*x1_6 = 1, b*x1_0*x4_1 = 1, -gama*sgm*x2_0*x4_0 = 1, x4_1 = 1, x3_0*x4_6 = 1, delta*sgm*x3_0*x4_7 = 1, -15*gama*sgm*x2_2*x4_4 = 1, -alpha*x1_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_3*c = 1, -gama*x4_0*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x3_0*x4_7 = 1, -7*gama*sgm*x2_1*x4_6 = 1, delta*sgm*x3_2*x4_4 = 1, x1_5*x4_4 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, -alpha*x1_5 = 1, b*x1_1*x4_3 = 1, x2_3 = 1, x3_2*x4_2 = 1, x1_5*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, x1_2*x4_1 = 1, delta*sgm*x3_1*x4_3 = 1, -gama*sgm*x2_0*x4_2 = 1, -gama*sgm*x2_0*x4_5 = 1, -x1_3 = 1, x1_4*x4_4 = 1, b*x1_1*x4_4 = 1, x3_5*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, b*x1_3*x4_4 = 1, x1_5*x4_1 = 1, b*c*x1_3 = 1, -gama*sgm*x2_0*x4_7 = 1, -gama = 1, -10*gama*sgm*x2_2*x4_3 = 1, -gama*sgm*x2_2*x4_0 = 1, x1_8*x4_0 = 1, b*x1_2*x4_5 = 1, x3_1*x4_6 = 1, -alpha*x1_6 = 1, b*x1_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, -x1_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_6*x4_3 = 1, x1_2*x4_3 = 1, b*x1_0*x4_0 = 1, -21*gama*sgm*x2_2*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, -alpha*x1_2 = 1, x1_2*x4_0 = 1, -gama*sgm*x2_0*x4_3 = 1, x1_4*c = 1, x2_5 = 1, -gama*x2_2 = 1, -x1_1 = 1, x4_2 = 1, -1781155436326284518167226492020936932495795105789247131042693581414951903774640915090324086302795368207309518909758946768057268927599549017003023092708938652628413809138/415360464679697589169918132066320839991110185956005680411835121473290696456549495792785198165 = 1, -x1_7 = 1, delta*sgm*x3_1*x4_6 = 1, x1_7*c = 1, b*x1_2*x4_4 = 1, -2*gama*sgm*x2_1*x4_1 = 1, x1_3*x4_5 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_1*x4_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x3_1*x4_7 = 1, x1_6*x4_2 = 1, b*x1_1*x4_0 = 1, b*x1_3*x4_5 = 1]
# 282, -5.588688399
# 2
# [x3_4 = 10, x2_4 = 7]
# [x3_4 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2], x2_4 = [4, 1, 4, 2, 2, 4, 4]]