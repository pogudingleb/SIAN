infolevel[Groebner]:=10;
Et_hat := [403095681992776830-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-22474048242896632629078443436738225858501681779352140789/12166328170900227333, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+935566807646739260962304188251422443297180798250469151818825239758537242255127124236710029973909958162552577579/110519621113067482652040192126125465694499316487979750395, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-7789297612164395282175786699358411192638597159123757126233612370496884161166506359069403955077738650637953430221532698666110269140452173803284547117427653221168138193/200793312154626740330864984123962859057992977813953448009155894648230850031727345459492976385, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+45292303832699051997007618506177683429380134216099321165096194802916230340938221185738022732063535510328261039304528466087079101001205514783015313586234721274920966000885078008413330042787341924898477875206209886903910981/608005978183999735502658998938851680410145081472822706195317907245612386876668631655395425241309120910612223556437722861415427925, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5-10714054045396434362418415841739198287313949590382083623663522353865055613818184065563277087239509899841986893825220324605306299941052541790783025104423145540783677246191474358904469268569398444079150383345384976823299852322922103363323409173785683448494478124195017973543247961621974547841031483/5523161088494913075000925450481751317626962262777564681984927862265367562882063348692190993736558762944330427487953867691093069081864713303886117521557851308514408875, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+2091426233437974680892100726317882607425946459297164867415546437974143632776020125254469986991892144445354671375672501224729786586951661375316369958335467037486963395536344759230734234514175787634320761948885163796197638133631181122635415454515589233983073202347368654245612938479831373802799219018459366108862538881052773560450132159409343476457423128760290621645803621/50172711295665167640779924449073559291341333754083891719973712534978018639636211690063566312379108069324239100658070062872937069613816319175354920571327016370282603742720634957297788618108097711572983125, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7+105774762804938322899825168371365048003645468273611849529424832188134858681844250341177666559181449055829363875868171908339545203632515660289276988492461964621702911668872028356262715435273674838715722088770233327848004454962449895266175352676534483371124375724936621472086961771847591277081675707616502870409496442126542908467625647357213150539775547988246027599814145484609686965839247666262347637227402934137503740758163756789072291383994301/50641308853243099181881364761754532028923669240717342221585324060140929074035482896665910306358847835488345603382907001698003673327987919053889177737152368409050259925548120357179262005717796751380879844930716702600061710040545608469721875, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-6938149775993575374245454774250022418258354304216234958764222925750755255467525106725857885297173748304124148203578362849840782666546775789308354330158638549241364393788975230156633172676666259681463591327007073342895446914831943032934115152165484275163975962773290483804091915022770533232574879470298440055947692833500203204784737957163453372554307400899005689657221005931081803334180421682419760593467155796777549434970185203050075031029881078393188760253786885527042807270170190714453386818561263055060766331947093/276017128019790722823564230239625792998282292638361834805231388291871421297475054973045566964021292620758716798784195418405203407929815417308072106179610467404705499731403603511753785646158553273447463171253010026492401218329902920139684130144545416246501824886817401240246875, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9-9254581255355786418859588547198131215942755514134686943486160924968772967132964776516268493005806046626359708827069610788676805186932945951450098919273192151821450455094414240833877822039203251283331858660628466325453811774808310978513067139862944250963859904261564332568694829705562076898594072836824625516639398661262363280971137562242864736482867847155240521062650897421981690942432861557625412577045957665114912063423394143162361369650027588810030945676354342718526085160756640729165054126559075349454555999243502098183093731439425426408487345501407562776552184117540877395626972029239381/2507355381258561199246889971419447726091702628837211572388776827624116527788232598301765683888785354459066646077854520179823244817309983042283559834822329715302993331678012798922675027469345889232015308016400334333816727646500273978326269963186672572768227225516498722325122154928849077332322136566790273494953125, 9476545671845013374967341650444914637548663712364048358031633994167042711262756980665647723474151779846504007665724053166285931058/110519621113067482652040192126125465694499316487979750395-x2_3, 50658446236408232561188790241433842277175900142234443456-x3_2, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_3, x3_2] 64
# 274, [1 = 8, -1 = 2, delta*sgm*x3_0*x4_3 = 1, b*x1_6*x4_0 = 1, x4_1 = 1, x3_7 = 1, x3_3*x4_1 = 1, x4_5*delta*sgm = 1, x1_5*x4_1 = 1, b*c*x1_1 = 1, -15*gama*sgm*x2_2*x4_4 = 1, x1_5*c = 1, x3_0*x4_5 = 1, b*x1_2*x4_1 = 1, x1_7*x4_2 = 1, z_aux*x3_0*c = 1, x3_6*x4_1 = 1, x3_3 = 1, b*x1_2*x4_2 = 1, x4_6 = 1, -gama*x2_4 = 1, delta*x3_4 = 1, x2_6 = 1, b*x1_2*x4_5 = 1, -gama*sgm*x2_0*x4_1 = 1, x3_0*x4_1 = 1, b*x1_5*x4_3 = 1, -20*gama*x4_3*sgm = 1, delta*sgm*x3_6*x4_0 = 1, b*x1_4*x4_0 = 1, x1_5*x4_3 = 1, x1_1*x4_2 = 1, x1_1*x4_8 = 1, delta*sgm*x3_4*x4_0 = 1, x1_6*x4_0 = 1, x2_7 = 1, -x4_0*gama*sgm = 1, delta*sgm*x3_4*x4_3 = 1, z_aux*x3_0*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, -x1_7 = 1, delta*x3_3 = 1, x1_8*x4_1 = 1, b*x1_5*x4_1 = 1, -x1_6 = 1, -6*gama*sgm*x2_1*x4_5 = 1, b*x1_0*x4_2 = 1, -10*gama*x4_2*sgm = 1, x3_1*x4_4 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_3*x4_3 = 1, b*x1_5*x4_0 = 1, x4_4*delta*sgm = 1, b*x1_1*x4_2 = 1, -3*gama*sgm*x2_2*x4_1 = 1, x2_2 = 1, x1_6*x4_2 = 1, x1_1*x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, -gama*sgm*x2_0*x4_4 = 1, x1_4*c = 1, x2_4 = 1, x4_3*delta*sgm = 1, x3_1 = 1, delta*sgm*x3_1*x4_6 = 1, x3_1*x4_6 = 1, x3_1*x4_5 = 1, x1_4*x4_3 = 1, x3_0*x4_4 = 1, -x1_2 = 1, x3_3*x4_4 = 1, delta*sgm*x3_0*x4_1 = 1, x1_3*x4_5 = 1, -3*gama*sgm*x2_1*x4_2 = 1, delta*x3_0 = 1, x3_6*x4_2 = 1, -6*gama*sgm*x2_5*x4_1 = 1, delta*sgm*x3_1*x4_4 = 1, x1_1*c = 1, x1_2*x4_1 = 1, b*x1_1*x4_6 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_1 = 1, b*x1_0*x4_5 = 1, x1_6*c = 1, x1_5*x4_4 = 1, beta*x2_2 = 1, x3_3*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x1_4*x4_0 = 1, -gama*sgm*x2_6*x4_0 = 1, x3_4*x4_2 = 1, x4_2 = 1, -35*gama*x4_4*sgm = 1, delta*sgm*x3_5*x4_1 = 1, x3_4*x4_3 = 1, delta*sgm*x3_0*x4_4 = 1, x1_2*x4_0 = 1, b*x1_7*x4_1 = 1, b*x1_4*x4_3 = 1, delta*x3_5 = 1, x1_9*c = 1, -gama*sgm*x2_2*x4_0 = 1, -15*gama*sgm*x2_4*x4_2 = 1, b*c*x1_7 = 1, delta*sgm*x3_3*x4_1 = 1, -alpha*x1_3 = 1, b*x1_0*x4_1 = 1, x1_8*x4_0 = 1, x1_4*x4_2 = 1, delta*sgm*x3_3*x4_4 = 1, x1_5*x4_0 = 1, x1_1*x4_4 = 1, x1_1*x4_3 = 1, x3_1*x4_1 = 1, x1_1*x4_0 = 1, beta*x2_6 = 1, -alpha*x1_4 = 1, x3_5*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, beta*x2_0 = 1, x2_5 = 1, -gama*sgm*x2_0*x4_2 = 1, x3_5*x4_1 = 1, x1_5*x4_2 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_7 = 1, x2_1 = 1, -gama*sgm*x2_0*x4_0 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*c*x1_6 = 1, -x1_4 = 1, b*c*x1_3 = 1, x3_0*x4_8 = 1, x3_7*x4_1 = 1, b*x1_4*x4_2 = 1, -5242672485478177444659387929538951251864037998431138078655813910570763801400630223669703878942695805601358660228072931893481643120033360259685914626310457154998773153/3446941865585422159032084916215751002267832245135820694126174656986797068950014823185448279 = 1, x1_6*x4_3 = 1, b*x1_7*x4_0 = 1, x3_0*x4_2 = 1, -7*gama*sgm*x2_6*x4_1 = 1, b*x1_2*x4_3 = 1, b*x1_5*x4_2 = 1, delta*sgm*x3_0*x4_5 = 1, x1_2*x4_3 = 1, x1_6*x4_1 = 1, -gama*sgm*x2_4*x4_0 = 1, x1_3*x4_4 = 1, b*x1_1*x4_7 = 1, delta*sgm*x3_4*x4_1 = 1, b*c*x1_0 = 1, x4_3 = 1, x4_5 = 1, delta*sgm*x3_0*x4_2 = 1, -gama*x2_1 = 1, x1_2*c = 1, -gama*sgm*x2_0*x4_7 = 1, x1_2*x4_4 = 1, b*x1_1*x4_1 = 1, -15115441682291351469493999796889194256444309189702658997423903023753904476681805540674372864157150361997078285344947881840603112801189118660874555200030360184490623136736094978269849325660690618500256904101614451942602254596539800869941050913585838514524194595778754542279355296272339717990213237690915121173658700936742021099515999046296336730983160044154731834/55966407153616880985053572623591423628062510465785284303202635187864245602512659742542807287088833506264938261188496067030530498949786507070432964268747772403978351984641552701150465895475535612767 = 1, x1_2*x4_2 = 1, x3_5 = 1, -alpha*x1_1 = 1, -285178703914060845731182048885212093993751879847955314363/8614621920582149333 = 1, x4_2*delta*sgm = 1, beta*x2_1 = 1, b*x1_0*x4_7 = 1, x1_3*x4_0 = 1, x3_4 = 1, delta*sgm*x3_0*x4_0 = 1, x1_2*x4_6 = 1, delta*sgm*x3_6*x4_1 = 1, x1_8*c = 1, x1_3*x4_3 = 1, delta*x3_6 = 1, x3_4*x4_1 = 1, b*c*x1_4 = 1, x1_9*x4_0 = 1, -gama = 1, beta = 1, x1_3*x4_2 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_0*x4_6 = 1, beta*x2_5 = 1, -alpha*x1_2 = 1, x3_5*x4_3 = 1, delta*sgm*x3_4*x4_2 = 1, b*x1_0*x4_6 = 1, delta*sgm*x3_3*x4_0 = 1, x1_4*x4_5 = 1, b*x1_4*x4_1 = 1, -x1_8 = 1, delta*x3_1 = 1, b*x1_3*x4_0 = 1, -4*gama*x4_1*sgm = 1, -gama*x2_0 = 1, -x1_9 = 1, b*c*x1_8 = 1, b*x1_0*x4_8 = 1, b*x1_6*x4_1 = 1, -5629317267317686165489021192098964322447737016015043776670822618707790786426718509586239677710877917116627651791041876428312678539839048764250576436518854083802900965418286907046686072707809546263116255395066494176408946276651216573074061161588942758327385018249058296052529508589053836325/153246013764087925164908048586910532790333032454155717482009559419018673359698329123527617565138776176515502005311383925632499363141958349866424164716281220390053 = 1, x1_3*x4_1 = 1, -4*gama*sgm*x2_1*x4_3 = 1, delta*sgm*x3_3*x4_2 = 1, -gama*x2_6 = 1, x3_6 = 1, x1_1*x4_7 = 1, x1_7*x4_0 = 1, x3_1*x4_2 = 1, -x1_3 = 1, -x1_5 = 1, -gama*x2_2 = 1, delta*sgm*x3_1*x4_1 = 1, x3_3*x4_5 = 1, -10*gama*sgm*x2_2*x4_3 = 1, b*x1_2*x4_4 = 1, b*x1_4*x4_4 = 1, delta*sgm*x3_1*x4_2 = 1, -21*gama*sgm*x2_2*x4_5 = 1, x1_3*x4_6 = 1, b*x1_0*x4_0 = 1, b*x1_6*x4_2 = 1, x1_1*x4_6 = 1, -35*gama*sgm*x2_4*x4_3 = 1, x3_0*x4_3 = 1, x1_4*x4_1 = 1, -x1_0 = 1, -alpha*x1_0 = 1, delta*sgm*x3_5*x4_2 = 1, -gama*sgm*x2_1*x4_0 = 1, delta*x4_1*sgm = 1, b*x1_1*x4_3 = 1, b*x1_1*x4_0 = 1, x3_0*x4_7 = 1, x4_4 = 1, -5*gama*sgm*x2_1*x4_4 = 1, x3_4*x4_4 = 1, delta*x4_0*sgm = 1, x3_1*x4_7 = 1, b*x1_8*x4_0 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*x1_0*x4_4 = 1, x1_2*x4_5 = 1, -alpha*x1_6 = 1, x1_2*x4_7 = 1, x1_7*c = 1, b*x1_1*x4_4 = 1, x3_1*x4_3 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*x1_3*x4_4 = 1, x1_3*c = 1, -gama*x2_5 = 1, b*c*x1_5 = 1, b*x1_3*x4_5 = 1, delta*sgm*x3_3*x4_3 = 1, b*x1_1*x4_5 = 1, b*x1_0*x4_3 = 1, -alpha*x1_5 = 1, beta*x2_4 = 1, -21*gama*sgm*x2_5*x4_2 = 1, b*x1_2*x4_6 = 1, x1_7*x4_1 = 1, -gama*sgm*x2_0*x4_5 = 1, x1_1*x4_5 = 1, delta = 1, b*x1_2*x4_0 = 1, b*x1_3*x4_2 = 1, -x1_1 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_3*x4_3 = 1, -6*gama*sgm*x2_2*x4_2 = 1, delta*sgm*x3_0*x4_6 = 1, b*c*x1_2 = 1, -gama*sgm*x2_5*x4_0 = 1]
# 282, -5.577999884
# 2
# [x2_3 = 8, x3_2 = 14]
# [x2_3 = [4, 1, 4, 2, 2, 4, 4, 4], x3_2 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2]]
# [x2_3 = [1, 1, 1, 1, 1, 1, 1, 1], x3_2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]