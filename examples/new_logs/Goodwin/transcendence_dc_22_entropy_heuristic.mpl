infolevel[Groebner]:=10;
Et_hat := [1565953734535765216-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, 62054105181788637-gama_0, gama_1, -x1_1-129076934724629043956922161076690119869305961292387371103/14242135200075420129, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama_0*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+72383915659971186799083638310403724633016970757785720617149878580742116519031202490734249251490507961740365206149/1379983587883802751730484146917638532411894911802227874225, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_1-gama_0*x2_1-gama_1*x2_0)*x4_0-x4_1*(-delta*x3_0+gama_0*x2_0))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama_0*x2_0+x3_1, -x1_3-40591537577550571068202116457252783427921419752114231866323898965300916575097682912111527664336766714574200191873220197286765181869600962051359565999809507350526535717463/133712724677586863653982877505682057407943758625851444072035635765301018964337748820209385130625, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0)*x4_0+(2*delta*x3_1-2*gama_0*x2_1-2*gama_1*x2_0)*x4_1+x4_2*(delta*x3_0-gama_0*x2_0))*sgm+2*x3_1*x4_2+x3_0*x4_3+x3_2*x4_1, gama_2, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama_0*x2_1-gama_1*x2_0+x3_2, -x1_4+7570357579136131438368619822957528225398629068772303208772379507534990564059700396013165023635012560490584532459860080819615752915125206221602489533855596047024703204297124704640227602681035160005728685653305083749839623608193/4318672796227873379607064609310546339015469589386527000372296178917751431988897601235531472500251281872074935262089507501692020046875, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0)*x4_0+(3*delta*x3_2-3*gama_0*x2_2-6*gama_1*x2_1-3*gama_2*x2_0)*x4_1+(-x2_0*x4_3-3*x2_1*x4_2)*gama_0+(x3_0*x4_3+3*x3_1*x4_2)*delta-3*x4_2*gama_1*x2_0)*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, gama_3, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama_0*x2_2-2*gama_1*x2_1-gama_2*x2_0+x3_3, -x1_5+434350489731202846424261913721101190655041601522017758726666416818180655792326774805749880797221860087888828362470275041452059696554258061981734591891414160281641934846466519370847034090783373725970297235214731293932233023125662760849367448469480238383830284608828332993828264276315450451529405691/18193710595851482941231554783797771122092770672389907438478184667440988429766460375202578946287578440324746262915259078761319923780044699935611871180020318405234680078125, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0)*x4_0+(4*delta*x3_3-4*gama_0*x2_3-12*gama_1*x2_2-12*gama_2*x2_1-4*gama_3*x2_0)*x4_1+(-x2_0*x4_4-4*x2_1*x4_3-6*x2_2*x4_2)*gama_0+(x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2)*delta+(-4*gama_1*x4_3-6*gama_2*x4_2)*x2_0-12*x4_2*gama_1*x2_1)*sgm+4*x3_3*x4_2+x3_4*x4_1+6*x3_2*x4_3+4*x3_1*x4_4+x3_0*x4_5, gama_4, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama_0*x2_3-3*gama_1*x2_2-3*gama_2*x2_1-gama_3*x2_0+x3_4, -x1_6+2475798962240308576214583601372656595852528249077996992202512506674227466689760068473856031192753918486681767532302436110655672588753076385495702159104326962957321329893524488722466385901275722723865837009829982274371851346653392458084178138111163101388995739703813685522115980372762168993696603045462694396837113620253727069574708192992925677040867964598648589928851244039/8109198493938519070871766940601886371768450314839011367865362477009878701362191436929145674194671263713547160413337594131458374807677827768539623902925962327134365139847999300658111733490815204443358287109375, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0)*x4_0+(5*delta*x3_4-5*gama_0*x2_4-20*gama_1*x2_3-30*gama_2*x2_2-20*gama_3*x2_1-5*gama_4*x2_0)*x4_1+(-x2_0*x4_5-5*x2_1*x4_4-10*x2_2*x4_3-10*x2_3*x4_2)*gama_0+(x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2)*delta+(-5*gama_1*x4_4-10*gama_2*x4_3-10*gama_3*x4_2)*x2_0+(-30*gama_1*x2_2-30*gama_2*x2_1)*x4_2-20*x4_3*gama_1*x2_1)*sgm+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+5*x3_1*x4_5+x3_0*x4_6, gama_5, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama_0*x2_4-4*gama_1*x2_3-6*gama_2*x2_2-4*gama_3*x2_1-gama_4*x2_0+x3_5, -x1_7-230716284612356211720569419089056277317509866595030506166328918807499701969376049626640374652708654268476586988931088544093105285018625499715249115218489252466187943748235401226285626426879698369011888112230027682277644753814933247955402444171242114854534264873228246471423968076228164771164278142139638095532830056544728545570244900426987944896687476856200753818403152123980129265754811891932256777950625509082694484233533033388272953817587905479/52382411650691004648894199542993286409591848444498884722642924439251024860493350568223774524974105988400144532501374966929557913847581405603438607571550507046657570729315834057400810322248032310023773890176946091566514465635937152729562900390625, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0)*x4_0+(6*delta*x3_5-6*gama_0*x2_5-30*gama_1*x2_4-60*gama_2*x2_3-60*gama_3*x2_2-30*gama_4*x2_1-6*gama_5*x2_0)*x4_1+(-x2_0*x4_6-6*x2_1*x4_5-15*x2_2*x4_4-20*x2_3*x4_3-15*x2_4*x4_2)*gama_0+(x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2)*delta+(-6*gama_1*x4_5-15*gama_2*x4_4-20*gama_3*x4_3-15*gama_4*x4_2)*x2_0+(-60*gama_1*x2_3-90*gama_2*x2_2-60*gama_3*x2_1)*x4_2+(-30*x2_1*x4_4-60*x2_2*x4_3)*gama_1-60*gama_2*x2_1*x4_3)*sgm+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1+x3_0*x4_7, gama_6, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama_0*x2_5-5*gama_1*x2_4-10*gama_2*x2_3-10*gama_3*x2_2-5*gama_4*x2_1-gama_5*x2_0+x3_6, -x1_8-156137404851787483183307199664120131167959693921051743387361344441162754426535864265768607207755501277140816461614855561004842366124026108320607224891664131087260593085380490708343638256604142774711185070272640164903543968128013702446572942323026674763783228587421743426335985668864257426472612045793540966677943287801071068842621454263203366913617781860803094769655735771478903030821487860266078554815394826765183745403757395198763803545245230949193752485803832243406515747217750240541406787611833001029161935523420815303/3625402840446488443736680644335580152230967959956265438840133160760571338334761309620703826570938292501440728926509119798996863381256659111356423164574669045580986696260042197080951543911204307027578128698319608233627217357312857177131395626456014994349884109245778720305106201171875, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((delta*x3_7-gama_0*x2_7-7*gama_1*x2_6-21*gama_2*x2_5-35*gama_3*x2_4-35*gama_4*x2_3-21*gama_5*x2_2-7*gama_6*x2_1-gama_7*x2_0)*x4_0+(7*delta*x3_6-7*gama_0*x2_6-42*gama_1*x2_5-105*gama_2*x2_4-140*gama_3*x2_3-105*gama_4*x2_2-42*gama_5*x2_1-7*gama_6*x2_0)*x4_1+(-x2_0*x4_7-7*x2_1*x4_6-21*x2_2*x4_5-35*x2_3*x4_4-35*x2_4*x4_3-21*x2_5*x4_2)*gama_0+(x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2)*delta+(-7*gama_1*x4_6-21*gama_2*x4_5-35*gama_3*x4_4-35*gama_4*x4_3-21*gama_5*x4_2)*x2_0+(-105*gama_1*x2_4-210*gama_2*x2_3-210*gama_3*x2_2-105*gama_4*x2_1)*x4_2+(-42*x2_1*x4_5-105*x2_2*x4_4-140*x2_3*x4_3)*gama_1+(-105*gama_2*x4_4-140*gama_3*x4_3)*x2_1-210*gama_2*x2_2*x4_3)*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, gama_7, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama_0*x2_6-6*gama_1*x2_5-15*gama_2*x2_4-20*gama_3*x2_3-15*gama_4*x2_2-6*gama_5*x2_1-gama_6*x2_0+x3_7, -x1_9+476965809814831108050195216881216076414355941195096425264534775513729411740895570920642274745495764825943003768157605551767544576152660563274166810776167941087472714117648609789207704340169436505739941930354558721565521164971291574299839516412793771468892466390596512680179110398148042803223593182084414598319148913890315968772443527041638105566632056072410582759434505054488435037394095855377459497607408592847988259100050446702424709857077471086238399303687986466110560419174202406392492828247053195584908742503267153343188469036058738701620960153233954255978495233150501583555795433837753914653/351281345739294844374807007695898438447241811271220713755805980540103765398402864527077939508968894168599313603518271912979119239769684989852878858525567881078381049591461940102881476225718590806952123139362819181751875481016804087795857337763188636648223892185885348359880994681341080256851301550784933905524383544921875, -gama_1, -gama_2, -gama_3, -gama_4, -gama_5, -gama_6, -gama_7, -250821778683485959179784210435272274535272436099232039469041853057239417840297968757015530464274860332116842190309698037818289488358402596199223648148677561024679108634785501779324952555406/133712724677586863653982877505682057407943758625851444072035635765301018964337748820209385130625-x2_4, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, gama_7, x4_7, x3_7, x2_7, x1_7, gama_6, x4_6, x3_6, x2_6, x1_6, gama_5, x4_5, x3_5, x2_5, x1_5, gama_4, x4_4, x3_4, x2_4, x1_4, gama_3, x4_3, x3_3, x2_3, x1_3, gama_2, x4_2, x3_2, x2_2, x1_2, gama_1, x4_1, x3_1, x2_1, x1_1, gama_0, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [gama, x2_4] 179
# 274, [1 = 7, -1 = 3, x3_1*x4_5 = 1, delta*x3_2 = 1, x2_2 = 1, x3_1 = 1, x2_5 = 1, x1_3*x4_1 = 1, -21*sgm*x2_5*x4_2 = 1, b*x1_0*x4_5 = 1, -4*sgm*x2_3*x4_1 = 1, x1_2*x4_7 = 1, x3_4*x4_3 = 1, x1_1*x4_0 = 1, x1_1*x4_5 = 1, delta*x3_0 = 1, delta*sgm*x3_0*x4_1 = 1, b*x1_5*x4_0 = 1, x1_7*c = 1, delta*sgm*x3_0*x4_0 = 1, -x1_5 = 1, b*x1_2*x4_3 = 1, x2_7 = 1, x1_2*x4_5 = 1, -sgm*x2_5*x4_0 = 1, x2_6 = 1, -x1_9 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, delta*sgm*x3_2*x4_4 = 1, x3_0*x4_3 = 1, z_aux*x3_0*x4_0 = 1, x3_4 = 1, -15*sgm*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x3_1*x4_4 = 1, b*c*x1_7 = 1, b*x1_1*x4_3 = 1, x3_5*x4_2 = 1, delta*x3_6 = 1, x1_3*x4_6 = 1, b*x1_6*x4_0 = 1, delta*sgm*x3_3*x4_0 = 1, x3_0*x4_5 = 1, b*x1_6*x4_2 = 1, b*x1_7*x4_0 = 1, b*x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*x4_5 = 1, -alpha*x1_2 = 1, -x2_3 = 1, x1_6*x4_0 = 1, delta*sgm*x3_1*x4_4 = 1, -x1_4 = 1, x1_1*c = 1, x1_1*x4_6 = 1, -sgm*x2_2*x4_0 = 1, delta*sgm*x3_0*x4_3 = 1, beta*x2_3 = 1, x1_1*x4_1 = 1, -sgm*x2_0*x4_1 = 1, b*x1_3*x4_2 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, x3_7 = 1, delta*sgm*x3_6*x4_0 = 1, b*c*x1_3 = 1, x1_8*c = 1, -x1_2 = 1, b*x1_0*x4_3 = 1, delta*sgm*x3_4*x4_3 = 1, x3_1*x4_7 = 1, delta*sgm*x3_0*x4_7 = 1, x1_1*x4_7 = 1, x3_0*x4_7 = 1, x1_7*x4_1 = 1, x1_2*x4_3 = 1, x1_4*x4_3 = 1, x1_4*x4_5 = 1, -x2_2 = 1, b*x1_3*x4_0 = 1, beta*x2_5 = 1, -10*sgm*x2_2*x4_3 = 1, delta*sgm*x3_3*x4_2 = 1, -sgm*x2_0*x4_6 = 1, x1_2*x4_4 = 1, x1_1*x4_3 = 1, delta*x3_1 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_0*x4_6 = 1, x3_3*x4_2 = 1, -3*sgm*x2_1*x4_2 = 1, x1_3*x4_4 = 1, -sgm*x2_0*x4_3 = 1, delta*sgm*x3_1*x4_3 = 1, x1_1*x4_4 = 1, -35*sgm*x4_3 = 1, -x1_0 = 1, delta*sgm*x3_2*x4_3 = 1, x1_2*x4_2 = 1, x2_3 = 1, x1_9*x4_0 = 1, b*x1_7*x4_1 = 1, x3_6 = 1, b*x1_4*x4_3 = 1, -x1_6 = 1, b*x1_3*x4_5 = 1, x2_1 = 1, delta*sgm*x3_4*x4_0 = 1, x1_2*x4_0 = 1, b*x1_0*x4_1 = 1, b*x1_1*x4_0 = 1, b*x1_2*x4_2 = 1, b*c*x1_2 = 1, x3_2*x4_6 = 1, x3_2*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, delta*sgm*x3_7*x4_0 = 1, b*c*x1_4 = 1, b*x1_2*x4_0 = 1, -x1_8 = 1, x3_6*x4_1 = 1, delta*sgm*x3_5*x4_2 = 1, x3_1*x4_3 = 1, delta*sgm*x3_1*x4_2 = 1, delta*x3_5 = 1, x1_5*x4_4 = 1, -7*sgm*x2_6*x4_1 = 1, x1_5*x4_3 = 1, -21*sgm*x2_2*x4_5 = 1, -4*sgm*x2_1*x4_3 = 1, -alpha*x1_5 = 1, delta*sgm*x3_5*x4_1 = 1, x3_3*x4_5 = 1, -sgm*x2_0*x4_5 = 1, b*c*x1_5 = 1, -6*sgm*x2_5*x4_1 = 1, -x1_7 = 1, -3*sgm*x2_2*x4_1 = 1, x1_5*x4_2 = 1, x3_1*x4_6 = 1, b*x1_0*x4_7 = 1, b*x1_1*x4_2 = 1, x3_0*x4_6 = 1, x1_3*c = 1, x1_4*x4_2 = 1, x1_5*x4_0 = 1, -15*sgm*x2_2*x4_4 = 1, b*x1_4*x4_2 = 1, x3_2*x4_2 = 1, x3_1*x4_2 = 1, x3_4*x4_2 = 1, x3_4*x4_1 = 1, x1_1*x4_2 = 1, -5*sgm*x4_1 = 1, x1_6*x4_3 = 1, x1_9*c = 1, x3_3 = 1, beta*x2_0 = 1, x3_0*x4_8 = 1, b*x1_3*x4_1 = 1, -alpha*x1_0 = 1, x1_8*x4_0 = 1, -6*sgm*x2_1*x4_5 = 1, x3_0*x4_1 = 1, z_aux*x3_0*c = 1, beta = 1, -sgm*x2_0*x4_4 = 1, -5*sgm*x2_1*x4_4 = 1, b*x1_5*x4_1 = 1, b*x1_1*x4_7 = 1, -x4_0*sgm = 1, -sgm*x2_0*x4_2 = 1, b*x1_4*x4_4 = 1, x1_2*x4_6 = 1, -sgm*x2_1*x4_0 = 1, b*x1_5*x4_2 = 1, x1_3*x4_3 = 1, b*x1_0*x4_8 = 1, delta*sgm*x3_3*x4_1 = 1, x1_3*x4_2 = 1, -x2_0 = 1, x1_3*x4_0 = 1, -sgm*x2_7*x4_0 = 1, -6*sgm*x2_2*x4_2 = 1, b*x1_1*x4_6 = 1, beta*x2_1 = 1, -sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_4 = 1, b*x1_1*x4_1 = 1, x3_5*x4_1 = 1, b*x1_1*x4_4 = 1, x3_1*x4_1 = 1, x3_3*x4_1 = 1, -7*sgm*x2_1*x4_6 = 1, -sgm*x2_6*x4_0 = 1, delta*x3_4 = 1, x1_2*x4_1 = 1, delta*sgm*x3_0*x4_6 = 1, x1_5*c = 1, x1_6*x4_2 = 1, b*x1_1*x4_5 = 1, delta*sgm*x3_1*x4_6 = 1, b*c*x1_1 = 1, x3_5 = 1, b*c*x1_6 = 1, x3_2*x4_5 = 1, x3_2*x4_4 = 1, delta*sgm*x3_1*x4_0 = 1, x1_8*x4_1 = 1, delta*sgm*x3_4*x4_2 = 1, -alpha*x1_4 = 1, delta*x3_3 = 1, x1_4*x4_1 = 1, x1_4*c = 1, delta*sgm*x3_2*x4_5 = 1, x1_1*x4_8 = 1, -2297306029824391595558031449960236824838226366293027466419041173826341843762215900541735527866801767603403013545732635182453132747756655164101319318982108026958082169588469820070381493073453267729147075450325641089459981159602020791921907946395140784584785846618023368797774669361300450016085425238991483784644280864603264397422678950382311065320364735644717259791515081650296884698530523574109896526731969865352999328482087931660360244433160265535043429420781305291702276797676789855187147807992622446327750748035447408/25057096371042314787586416995394769810497531634114477285846862462606603393646412137465103463968037118599488732942595890093010277896028612630792511439350000488397211447982301841892582503247439994451130958922556404167066405809035641254183621142085584562898685671827461576409918515625 = 1, x3_2*x4_1 = 1, b*x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, -35*sgm*x2_3*x4_4 = 1, -alpha*x1_3 = 1, -2*sgm*x2_1*x4_1 = 1, x3_7*x4_1 = 1, b*x1_0*x4_4 = 1, x3_4*x4_4 = 1, delta*sgm*x3_2*x4_1 = 1, -x2_1 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_0 = 1, b*x1_3*x4_3 = 1, -x1_3 = 1, -alpha*x1_6 = 1, delta*sgm*x3_2*x4_0 = 1, b*x1_4*x4_0 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_5 = 1, x1_4*x4_0 = 1, -10*sgm*x2_3*x4_2 = 1, b*x1_8*x4_0 = 1, beta*x2_2 = 1, -alpha*x1_1 = 1, b*x1_0*x4_0 = 1, x1_6*x4_1 = 1, b*x1_6*x4_1 = 1, b*x1_5*x4_3 = 1, -628698479548236224015748601630649329907573300944618305839/16046202655528427505 = 1, b*c*x1_8 = 1, delta*sgm*x3_2*x4_2 = 1, -x2_0*x4_0*sgm = 1, x3_3*x4_4 = 1, b*x1_2*x4_5 = 1, -198426303481760771902847908118347140290786563780517916406904359248859766770598661528540187913895309318170427123345915897671959050144191463885043671133466030542288187273077216640906596212743105633232030201478398934824796560247623324096266415737986735940671145470163004991855033895730908647594846954957378739129216016464753863210182792951961640463129054406605340813613371143236908506279662955974472700784746523157528291815562761151435744700452128/67141237653384923587123056390913328422680243394455544664430931317767623284230837163376721676944495007067434272873916581429167276376476037769462872117056210207404360873122534886844340631068032316439885788183283524349276571290796459682408859375 = 1, x3_3*x4_3 = 1, -20*sgm*x2_3*x4_3 = 1, x1_6*c = 1, -x1_1 = 1, delta*sgm*x3_5*x4_0 = 1, x3_0*x4_4 = 1, -39807190786230863282904089587560109583540687071500568010313300634951599389674272256774589331678477686786623900529631896602061279563381105548082895728572324021930414961462/27591166081327090329279230743899663049639823818727880781844072092272973197930329551816798532625 = 1, -x2_6 = 1, b*x1_0*x4_2 = 1, x1_5*x4_1 = 1, -x2_5 = 1, x3_5*x4_3 = 1, x1_7*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, -sgm*x2_0*x4_7 = 1, x1_7*x4_0 = 1, b*x1_2*x4_4 = 1, x3_0*x4_2 = 1, x3_6*x4_2 = 1, x3_2 = 1]
# 282, -5.581916971
# 2
# [x2_4 = 7, gama = 43]
# [x2_4 = [4, 1, 4, 2, 2, 4, 4], gama = [4, 4, 4, 2, 4, 4, 4, 2, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 2]]