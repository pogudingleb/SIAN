infolevel[Groebner]:=10;
Et_hat := [202298872397594419-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-24421570549318759471046942602837765736183745669612537951/14145374495973723832, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+1235558470941740080135112607079464920011843542752300416123497104539720446487062797973657423531588612540967464665/83856838371984235320565091334834360975092278378245061936, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-31255253056574517813607639945163083404910615371895769207753548667044758700551463734516414724332627428521197066679766132642321887398115576361297404260743887633304551499/248560734243785344078587316935615781517926255316508121331981554621028157050540516218142868464, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+6337949751455540531088147297267541910908158373445257121015061881886782882874613983851897665777378673296686345236319420896888069095481598323705373696064315836783532642663467593151733122064573424941483265358442836346616644569/5894087094840971616292996708905358399275464368513840824264454961177773044001297143452366730728940377106161651637352182572438089088, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-115534883147213162908661769756355434473559180532793729736957188257466682268848967765916742402483889571265988183983749281513586257179534378341618001663957858910287968327425387984719196224135086662871627734949719050304929968256363241371675775756459999022829665597306560262702651300362830943354603/34941422653968001397559074978991973159685257291530477219384659232440381824638982889720588079183995308087043183831760150460493950759357645469916348629008410238336788224, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6-10377676144329263117130364862528403811101650387668338713623831874159178929755478813613564572160736173096260239301962190382526666837021447083362279917974969917812495468512876568999681766825345005402639691507409500895493526217516981496756866534761317932373898108024978012287941050857936742240000783540623901282311188506564630484129100067533727193182249797802065787531039/103570154074570002493072933364279247194889920799475607674276624045140227538594295070548767108945355897244407848021421488260872076965341744707996149437967989658792936570697588630003105091454440085931205376, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+810288041118176215836161709664868975059638538768037181496580500082003282362670429518705652426783810341186362957983148754674263151671573249754754497411840123595171070721032407498758257787639198758679346000119436668020455060613889861100454049151727396362392014103606679441389382131129629226885708028787986280329326459095315878176784833816648092830603008757206480113454981013785495779423656746867724625390316638208742577467221555282912889491585/2455945064689464235860944189355524058867140034194291967846060101494239077542083560724998266527227021680058399160288833984138034958052589707854043148749829999799471988632504312875527899331288506366313669378990278414330188481091731798910728192, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8+141847909476356671264005512171122105657463453141585304411737991738132867392964400052240623339099059433190804851738867028348121693125788418439020254210726356141054646848730831888226727020213680878434217629400483347958488102471434871871976699521025472466203628467711657266805844266075416981909733609673048048175949244227411646681875644885797227718030709089824202219985567554075103946346750786063802251222991770579616131525711386497689694225707602619468685730016200704887345865908289026092515955200842848060124644308169/14559373341353163008622774661708054319971581450784743418043713670021958257283921896876302794125437626015414598073702228987437438713839954317719252120937114634414725881625591305264304373560492967826592794381551228038587358158340403335438306632648523301914095535813486853345161216, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-2430226427702000459189681476579760349536404850491511918275074800125159654033803212095293617331241296852562322818182637658023433291252399853398762509726124903031190183825552862458655185218017362253706619055647539214820251865979803110250425090031755349406378244417422249718145782200420723429581793501557458429365448952888516688084954242433792061296897182904587857831520665517886275795490312854201079790062671989514888477503234661232052857994490131070659795000236614336026983804378420328150498318268404343079159678454548442706771024529581607706408502058352184862650195839170214943313906629399/43155556518872675579885773306546494082481973246526340301477047840873259744108165490790312357651370991395603683704519368973897693657858831147943842580980672143538785531669225097723931129925811130877690641416432610296321568664464034146965049886724416856890905904584219646551564076999376976324785845848058461104861184, -27491053402137598175145883663392217100-x2_1, -94529750908880409176693511486552293970361018993872295094-x3_2, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_1, x3_2] 72
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, b*x1_6*x4_0 = 1, x4_1 = 1, x3_7 = 1, x3_3*x4_1 = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, x3_6*x4_1 = 1, x2_3 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, x4_6 = 1, -gama*x2_4 = 1, delta*x3_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, -3*gama*x4_2*sgm = 1, delta*sgm*x3_6*x4_0 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -gama*sgm*x2_3*x4_0 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, b*x1_0*x4_2 = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, b*x1_1*x4_2 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x2_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, x2_4 = 1, -2*gama*x4_1*sgm = 1, -gama*x2_3 = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x3_4*x4_2 = 1, x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -gama*sgm*x2_2*x4_0 = 1, -15*gama*sgm*x2_4*x4_2 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, -7*gama*x4_6*sgm = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, x2_5 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, -gama*sgm*x2_0*x4_0 = 1, -20*gama*sgm*x2_3*x4_3 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, x3_7*x4_1 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, -6*gama*x4_5*sgm = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, -gama*sgm*x2_4*x4_0 = 1, x1_3*x4_4 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, x4_3 = 1, x4_5 = 1, delta*sgm*x3_0*x4_2 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, -4*gama*x4_3*sgm = 1, x1_2*x4_2 = 1, x3_5 = 1, -5*gama*x4_4*sgm = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, x4_4*delta*sgm = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, x4_5*delta*sgm = 1, x4_1*delta*sgm = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, x1_3*x4_2 = 1, beta = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, -35*gama*sgm*x2_3*x4_4 = 1, x1_3*x4_1 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_1*x4_2 = 1, -x1_3 = 1, -10*gama*sgm*x2_3*x4_2 = 1, x4_2*delta*sgm = 1, -x1_5 = 1, -gama*x2_2 = 1, delta*sgm*x3_1*x4_1 = 1, x4_3*delta*sgm = 1, x3_3*x4_5 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, -21*gama*sgm*x2_2*x4_5 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_4*x4_3 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, x4_4 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, -gama*x2_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_1*x4_5 = 1, delta = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1, -gama*sgm*x2_5*x4_0 = 1]
# 282, -5.577999884
# 2
# [x2_1 = 10, x3_2 = 14]
# [x2_1 = [4, 1, 4, 2, 2, 4, 4, 4, 4, 4], x3_2 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2]]
# [x2_1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], x3_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]