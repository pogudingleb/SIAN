infolevel[Groebner]:=10;
Et_hat := [3110543295284625987-x1_0, b*c*x1_0+b*x1_0*x4_0+c*x1_1+x1_1*x4_0-1, -x1_1-22689292036594061680277161361179442764707505444188845060/11628085870174397663, ((c+x4_0)*b+x4_1)*x1_1+b*x1_0*x4_1+(c+x4_0)*x1_2, delta*sgm*x3_0*x4_0-gama*sgm*x2_0*x4_0+x3_0*x4_1, -x1_2+446391793948916394346048894961883431632777556777518805680414526482702489918084477269094758393396799399148412/364692660743718126679496640255445156096347102025077811, ((c+x4_0)*x1_2+x1_0*x4_2+2*x1_1*x4_1)*b+2*x1_2*x4_1+x1_1*x4_2+(c+x4_0)*x1_3, ((delta*x3_0-gama*x2_0)*x4_1+x4_0*(delta*x3_1-gama*x2_1))*sgm+x3_0*x4_2+x3_1*x4_1, -alpha*x1_0+beta*x2_0+x2_1, delta*x3_0-gama*x2_0+x3_1, -x1_3-254688571513070302794454504022496896925177483293866912092164047984809604326515445489081976605139238261115881394644873875904302556603313683711096274956825844514392/331698734449731094276148723533255194572557684923226247020575829491178635942735057447790843, ((c+x4_0)*x1_3+3*x1_1*x4_2+x4_3*x1_0+3*x1_2*x4_1)*b+3*x1_3*x4_1+x1_1*x4_3+3*x1_2*x4_2+(c+x4_0)*x1_4, ((x3_0*x4_2+2*x3_1*x4_1+x3_2*x4_0)*delta-2*gama*(x4_1*x2_1+1/2*x4_2*x2_0+1/2*x2_2*x4_0))*sgm+x3_2*x4_1+2*x3_1*x4_2+x3_0*x4_3, -alpha*x1_1+beta*x2_1+x2_2, delta*x3_1-gama*x2_1+x3_2, -x1_4+1090705928908466513523888609372140581031711556632653510275947126623300683472743729739934745792861652210222210287756775435088698306845897016242417135464744877938981343269970895509136227717128688331774168578988861216303608/3921967205880573048223942494564328660816495741072366868539420210304601099798668368465940479039660113313363211444611834265024167, ((c+x4_0)*x1_4+4*x1_3*x4_1+6*x1_2*x4_2+4*x1_1*x4_3+x4_4*x1_0)*b+4*x1_4*x4_1+x1_1*x4_4+4*x1_2*x4_3+6*x1_3*x4_2+(c+x4_0)*x1_5, ((x3_0*x4_3+3*x3_1*x4_2+3*x3_2*x4_1+x3_3*x4_0)*delta-gama*(x2_0*x4_3+3*x2_1*x4_2+3*x2_2*x4_1+x2_3*x4_0))*sgm+x3_3*x4_1+3*x3_2*x4_2+3*x3_1*x4_3+x3_0*x4_4, -alpha*x1_2+beta*x2_2+x2_3, delta*x3_2-gama*x2_2+x3_3, -x1_5+97197950172777654092699953263928501075709918181469511144330228045768708691326280026032469640823481625942443323721307757569782700484505699124296536854283896702319759649887827957700681237691930393412742688910469849012712746700729428184861373148074253068850171003657086805877975815228212810293898/3567144883286062648791396194616724349988608802826516173617818615645091037783438130997130579331758904362574340824162415311511086556902753459748987751100277414856271, ((c+x4_0)*x1_5+5*x1_4*x4_1+10*x1_3*x4_2+10*x1_2*x4_3+5*x1_1*x4_4+x4_5*x1_0)*b+5*x1_5*x4_1+x1_1*x4_5+5*x1_2*x4_4+10*x1_3*x4_3+10*x1_4*x4_2+(c+x4_0)*x1_6, ((x3_0*x4_4+4*x3_1*x4_3+6*x3_2*x4_2+4*x3_3*x4_1+x3_4*x4_0)*delta-gama*(x2_0*x4_4+4*x2_1*x4_3+6*x2_2*x4_2+4*x2_3*x4_1+x2_4*x4_0))*sgm+6*x3_2*x4_3+4*x3_3*x4_2+x3_4*x4_1+4*x3_1*x4_4+x3_0*x4_5, -alpha*x1_3+beta*x2_3+x2_4, delta*x3_3-gama*x2_3+x3_4, -x1_6+551582105753086644266957874708242511492443944893743310341490550586261030934801083104070502305256518314330420678585747605327478575718008173669966437108097777015110476183647215533968642898043963131503481208730640562356256051245354949775598519550152016269592269724202145425915564647816323980948897229435185764133389529359904531267913955119827261559696345534187687431318/3244423512586966064286865094274916309136098445015696565396433675211539004591703334179560404183168048686124920926225190049140687806290135030855033381029391609375305156537836393114835081703311339408823, ((c+x4_0)*x1_6+6*x1_5*x4_1+15*x1_4*x4_2+20*x1_3*x4_3+15*x1_2*x4_4+6*x1_1*x4_5+x4_6*x1_0)*b+6*x1_6*x4_1+x1_1*x4_6+6*x1_2*x4_5+15*x1_3*x4_4+20*x1_4*x4_3+15*x1_5*x4_2+(c+x4_0)*x1_7, ((x3_0*x4_5+5*x3_1*x4_4+10*x3_2*x4_3+10*x3_3*x4_2+5*x3_4*x4_1+x3_5*x4_0)*delta-gama*(x2_0*x4_5+5*x2_1*x4_4+10*x2_2*x4_3+10*x2_3*x4_2+5*x2_4*x4_1+x2_5*x4_0))*sgm+5*x3_1*x4_5+10*x3_2*x4_4+10*x3_3*x4_3+5*x3_4*x4_2+x3_5*x4_1+x3_0*x4_6, -alpha*x1_4+beta*x2_4+x2_5, delta*x3_4-gama*x2_4+x3_5, -x1_7-3656885104265573882459864556526634560122324082858946741914156363858924351004015452404063844551136801975101527852269122906751468113215616313917608104276362185725011894204377943682680556774028192724073462361685192382403589417169754896631384823875566528073917544166586818367303951271885222343216366713976745675845322715346344559426334982741781924125016767950579131855794952150520655596169426909483116482346824404508158236071916100506270884077526/38361685761216968048841819250435839310421061159865717138437160931618195822993039695428908751706227319737538959120405803661276206878430797478555219691073047793947992367641602127565026037932249806076967154471568985328297107111841054824787, ((c+x4_0)*x1_7+7*x1_6*x4_1+21*x1_5*x4_2+35*x1_4*x4_3+35*x1_3*x4_4+21*x1_2*x4_5+7*x1_1*x4_6+x4_7*x1_0)*b+7*x1_7*x4_1+x1_1*x4_7+7*x1_2*x4_6+21*x1_3*x4_5+35*x1_4*x4_4+35*x1_5*x4_3+21*x1_6*x4_2+(c+x4_0)*x1_8, ((x3_0*x4_6+6*x3_1*x4_5+15*x3_2*x4_4+20*x3_3*x4_3+15*x3_4*x4_2+6*x3_5*x4_1+x3_6*x4_0)*delta-gama*(x2_0*x4_6+6*x2_1*x4_5+15*x2_2*x4_4+20*x2_3*x4_3+15*x2_4*x4_2+6*x2_5*x4_1+x2_6*x4_0))*sgm+x3_0*x4_7+6*x3_1*x4_6+15*x3_2*x4_5+20*x3_3*x4_4+15*x3_4*x4_3+6*x3_5*x4_2+x3_6*x4_1, -alpha*x1_5+beta*x2_5+x2_6, delta*x3_5-gama*x2_5+x3_6, -x1_8-6046011599983068408795761661967325284653641517885193807168353197227859142554630614468307085701727237378555671653890662945487828861784078559153290603832997703599857021657147344563339502410726556014027130911374123917779140210786128007806099458098031961277075535269106547678964598846283501923535167955395827807287280189352479751249480523127046336142935929270424339241314717097674914621059516014015679589584651930866407186103349020488080762382639591408339795140893215108631923780519381440081634808247782176626981241586/1203140861451412469743827528586771939509485869299398049207020831681381426421768789159631918878843806454520609111253477782764350665548008357200662023647317940258658870835209114924567928441624591777872376622820790817865631579158428113916883316585342463853432723195554674839, ((c+x4_0)*x1_8+8*x1_7*x4_1+28*x1_6*x4_2+56*x1_5*x4_3+70*x1_4*x4_4+56*x1_3*x4_5+28*x1_2*x4_6+8*x1_1*x4_7+x4_8*x1_0)*b+8*x1_8*x4_1+x1_1*x4_8+8*x1_2*x4_7+28*x1_3*x4_6+56*x1_4*x4_5+70*x1_5*x4_4+56*x1_6*x4_3+28*x1_7*x4_2+(c+x4_0)*x1_9, ((x3_0*x4_7+7*x3_1*x4_6+21*x3_2*x4_5+35*x3_3*x4_4+35*x3_4*x4_3+21*x3_5*x4_2+7*x3_6*x4_1+x3_7*x4_0)*delta-gama*(x2_0*x4_7+7*x2_1*x4_6+21*x2_2*x4_5+35*x2_3*x4_4+35*x2_4*x4_3+21*x2_5*x4_2+7*x2_6*x4_1+x2_7*x4_0))*sgm+x3_0*x4_8+7*x3_1*x4_7+21*x3_2*x4_6+35*x3_3*x4_5+35*x3_4*x4_4+21*x3_5*x4_3+7*x3_6*x4_2+x3_7*x4_1, -alpha*x1_6+beta*x2_6+x2_7, delta*x3_6-gama*x2_6+x3_7, -x1_9+456101752059367203353528280611803207251890185367431534571839770587959048407988631677086677789611444575830073811887098180986723144528451249380546725170463301761168702272432395411518734291036591790697817327497100738030744563152313942004182553071128640671266343502062501763268055930171794834455262490653728157361937580250878065000655289159503806838602094626836723625331136049177759885335149346658895701184679430314271510993195657727334294842085760544004279675455885605715462728796534120680992583790120550791763230372655726185635884664038357162235871952384342857160164732306069456737239984243/1094292109674891447135701382722918965189811013724197538081288231303858345225170389129070090521697781882205688863292967638834591696802863598266302133892689102138140624001622498395962645590800457382765021885926437517997608538502852395206180405933502642302256296894398860517730633798901590303419313219736188607, -53566789883231083824507336658530606756711680626022803882839882419071197241/11628085870174397663-x2_2, -28697055050413700094537475360888589847715268691117587269-x3_2, z_aux*x3_0*(c+x4_0)-1];
vars:=[x1_9, x4_8, x1_8, x4_7, x3_7, x2_7, x1_7, x4_6, x3_6, x2_6, x1_6, x4_5, x3_5, x2_5, x1_5, x4_4, x3_4, x2_4, x1_4, x4_3, x3_3, x2_3, x1_3, x4_2, x3_2, x2_2, x1_2, x4_1, x3_1, x2_1, x1_1, x4_0, x3_0, x2_0, x1_0, z_aux, w_aux, alpha, b, beta, c, delta, gama, sgm];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [x2_2, x3_2] 68
# 274, [1 = 8, -1 = 2, -7*gama*sgm*x2_6*x4_1 = 1, x1_2*x4_0 = 1, delta*x3_6 = 1, x1_2*x4_6 = 1, b*x1_3*x4_4 = 1, -x1_2 = 1, -15*gama*sgm*x2_4*x4_2 = 1, x3_6*x4_1 = 1, b*x1_1*x4_2 = 1, x3_1*x4_4 = 1, x3_1*x4_7 = 1, b*c*x1_2 = 1, x1_5*x4_4 = 1, b*x1_5*x4_2 = 1, x3_5*x4_2 = 1, x1_2*x4_4 = 1, -10*gama*x4_3*sgm = 1, delta*sgm*x3_1*x4_2 = 1, x1_2*x4_1 = 1, x1_1*x4_1 = 1, beta*x2_5 = 1, -gama*sgm*x2_0*x4_1 = 1, -1770977435351686016068346260438129682862046220309700163892609663766343204319395223811939121160074522582872295441873426109538757450368379806404807167743148460317914365946818078538984915278334890011404201915959662864877602272636609902307352381568554453912276230703458199280725004204607266229378349139204218727453582592021384244095365673331361848032317200583325478985536858394469179685595507363555448962462921753147816407357690229091755445755505/10293341692123280953159287302430384864146382227922150074030101361461568306292542099178095888330215093823260791239484669835274910692914088715490119593213342335562262684153258532947134992888627884734744720024482578963723021800411331065552 = 1, -x1_7 = 1, x1_5*x4_2 = 1, -20*gama*sgm*x2_3*x4_3 = 1, -gama*sgm*x2_5*x4_0 = 1, x1_8*x4_0 = 1, x1_1*x4_4 = 1, -gama*sgm*x2_0*x4_0 = 1, beta*x2_6 = 1, delta*sgm*x3_3*x4_4 = 1, x1_1*x4_3 = 1, b*x1_4*x4_2 = 1, -gama*sgm*x2_0*x4_2 = 1, b*x1_1*x4_3 = 1, delta*sgm*x3_0*x4_7 = 1, x3_0*x4_5 = 1, -2*gama*sgm*x2_1*x4_1 = 1, b*x1_4*x4_1 = 1, -gama*sgm*x2_1*x4_0 = 1, x4_4*delta*sgm = 1, delta*sgm*x3_0*x4_1 = 1, x1_1*x4_5 = 1, delta*sgm*x3_5*x4_0 = 1, b*x1_3*x4_2 = 1, x4_2*delta*sgm = 1, -gama*sgm*x2_3*x4_0 = 1, delta*sgm*x3_0*x4_0 = 1, -alpha*x1_6 = 1, x1_6*x4_3 = 1, delta*sgm*x3_1*x4_1 = 1, x1_9*x4_0 = 1, b*x1_0*x4_6 = 1, x1_7*c = 1, -21*gama*sgm*x2_5*x4_2 = 1, x4_3 = 1, -gama*x2_4 = 1, -alpha*x1_3 = 1, x4_2 = 1, x3_0*x4_3 = 1, x1_1*c = 1, b*x1_1*x4_0 = 1, delta*sgm*x3_7*x4_0 = 1, x3_4*x4_3 = 1, x1_3*x4_3 = 1, x1_7*x4_2 = 1, delta*x4_0*sgm = 1, x3_1*x4_1 = 1, delta*sgm*x3_4*x4_0 = 1, b*x1_2*x4_6 = 1, x1_4*x4_4 = 1, delta*sgm*x3_0*x4_4 = 1, -gama*sgm*x2_6*x4_0 = 1, b*x1_0*x4_2 = 1, b*c*x1_8 = 1, delta*x3_4 = 1, x1_7*x4_0 = 1, x1_3*x4_2 = 1, x1_2*x4_5 = 1, delta*sgm*x3_0*x4_3 = 1, x1_4*c = 1, b*x1_0*x4_3 = 1, b*x1_0*x4_8 = 1, x1_3*x4_4 = 1, delta*x3_1 = 1, b*x1_1*x4_6 = 1, -alpha*x1_4 = 1, -x1_5 = 1, x1_3*x4_0 = 1, x3_0*x4_8 = 1, x1_8*x4_1 = 1, b*x1_5*x4_3 = 1, beta*x2_1 = 1, -gama*sgm*x2_7*x4_0 = 1, x3_3*x4_4 = 1, b*x1_0*x4_1 = 1, beta*x2_4 = 1, delta*sgm*x3_4*x4_3 = 1, -gama*x2_3 = 1, x1_3*x4_5 = 1, x4_4 = 1, -3*gama*sgm*x2_1*x4_2 = 1, -272445264882798521927162215417260549350692976219005930595/7501903452127771988 = 1, x3_4 = 1, delta = 1, b*x1_1*x4_5 = 1, x1_4*x4_3 = 1, b*x1_5*x4_0 = 1, -6*gama*x4_2*sgm = 1, x1_8*c = 1, x1_2*x4_3 = 1, x4_3*delta*sgm = 1, delta*sgm*x3_6*x4_0 = 1, x1_2*c = 1, x1_1*x4_6 = 1, b*x1_2*x4_3 = 1, b*x1_3*x4_3 = 1, x1_1*x4_0 = 1, x3_3 = 1, delta*sgm*x3_0*x4_2 = 1, delta*sgm*x3_3*x4_3 = 1, -6*gama*sgm*x2_5*x4_1 = 1, b*x1_2*x4_2 = 1, b*c*x1_7 = 1, x1_5*x4_1 = 1, b*x1_6*x4_2 = 1, z_aux*x3_0*x4_0 = 1, x4_5 = 1, delta*sgm*x3_1*x4_3 = 1, b*x1_6*x4_0 = 1, x1_5*x4_3 = 1, b*x1_0*x4_0 = 1, delta*sgm*x3_5*x4_2 = 1, -7*gama*sgm*x2_1*x4_6 = 1, b*c*x1_1 = 1, delta*sgm*x3_1*x4_6 = 1, x1_5*c = 1, x1_3*x4_1 = 1, x1_2*x4_7 = 1, b*x1_7*x4_0 = 1, x3_6 = 1, -5*gama*sgm*x2_1*x4_4 = 1, beta*x2_0 = 1, b*x1_4*x4_3 = 1, x1_3*c = 1, b*x1_1*x4_1 = 1, x3_5*x4_3 = 1, delta*sgm*x3_3*x4_0 = 1, x1_1*x4_8 = 1, x2_4 = 1, x4_1 = 1, delta*sgm*x3_1*x4_0 = 1, x1_1*x4_2 = 1, -x1_8 = 1, x3_4*x4_4 = 1, x4_1*delta*sgm = 1, x1_4*x4_0 = 1, x3_5 = 1, b*x1_0*x4_7 = 1, delta*x3_3 = 1, -gama = 1, -gama*x2_1 = 1, -x1_1 = 1, -x1_6 = 1, -21*gama*x4_5*sgm = 1, -gama*sgm*x2_0*x4_7 = 1, delta*sgm*x3_0*x4_5 = 1, b*x1_3*x4_0 = 1, x3_1*x4_6 = 1, -gama*sgm*x2_0*x4_5 = 1, delta*x3_0 = 1, b*x1_2*x4_5 = 1, x1_5*x4_0 = 1, x1_3*x4_6 = 1, x3_0*x4_4 = 1, x3_1*x4_5 = 1, x3_3*x4_5 = 1, -3*gama*x4_1*sgm = 1, x3_6*x4_2 = 1, b*x1_1*x4_4 = 1, delta*sgm*x3_3*x4_2 = 1, -x1_0 = 1, delta*sgm*x3_4*x4_2 = 1, -gama*sgm*x2_0*x4_3 = 1, -x1_4 = 1, x2_1 = 1, b*x1_4*x4_0 = 1, x3_1*x4_2 = 1, delta*sgm*x3_1*x4_5 = 1, x3_3*x4_1 = 1, b*x1_2*x4_4 = 1, -gama*sgm*x4_0 = 1, x1_4*x4_2 = 1, -alpha*x1_1 = 1, b*x1_7*x4_1 = 1, -5*gama*sgm*x2_4*x4_1 = 1, b*c*x1_4 = 1, -4*gama*sgm*x2_1*x4_3 = 1, x2_5 = 1, x3_4*x4_2 = 1, b*x1_5*x4_1 = 1, -gama*x2_5 = 1, x3_0*x4_7 = 1, b*x1_8*x4_0 = 1, delta*sgm*x3_6*x4_1 = 1, b*c*x1_5 = 1, -35*gama*sgm*x2_4*x4_3 = 1, b*x1_1*x4_7 = 1, x1_1*x4_7 = 1, x1_6*c = 1, -gama*x2_6 = 1, b*c*x1_6 = 1, -4*gama*sgm*x2_3*x4_1 = 1, b*x1_3*x4_1 = 1, -x1_3 = 1, -35*gama*sgm*x2_3*x4_4 = 1, b*x1_6*x4_1 = 1, x1_6*x4_2 = 1, x3_1 = 1, x1_4*x4_5 = 1, -gama*sgm*x2_0*x4_6 = 1, x3_4*x4_1 = 1, x1_6*x4_1 = 1, b*x1_2*x4_0 = 1, b*x1_4*x4_4 = 1, beta = 1, -gama*x2_0 = 1, -alpha*x1_5 = 1, x3_0*x4_6 = 1, x1_4*x4_1 = 1, b*x1_2*x4_1 = 1, delta*x3_5 = 1, delta*sgm*x3_3*x4_1 = 1, x3_0*x4_2 = 1, delta*sgm*x3_4*x4_1 = 1, x4_5*delta*sgm = 1, x3_5*x4_1 = 1, x1_6*x4_0 = 1, x2_3 = 1, -alpha*x1_2 = 1, x3_7 = 1, x4_6 = 1, x3_7*x4_1 = 1, x3_1*x4_3 = 1, x3_3*x4_3 = 1, delta*sgm*x3_1*x4_4 = 1, x1_2*x4_2 = 1, -6*gama*sgm*x2_1*x4_5 = 1, -10*gama*sgm*x2_3*x4_2 = 1, -3336032149181473952864020087713184560238013976665987379171976632168321707062891453413701253431639255736689723519075308913578088832994438315617585162072130475330970534554109591586746796326537872800649354238067236193744601802894319479916635881922290489993337128746607006849023334969624596406736855010544450001221912630850208156451707132759842257636625339198715799064515/3075687755315850387175662323222007972494884287560494128503145195310479376820213705931964526464238216382270986848915206535415505804633662924467593798277116006413740347229059296594177936399730660633988 = 1, b*c*x1_3 = 1, -4556864563946229623878477120781517946238780291189388928946279646076171458135587639108897691164995222174383671048060933787619955927621795191956119414764057246930420885/2625728294855175635851655717123462957209652192301153216647477775470267752137313425729678394 = 1, -gama*sgm*x2_0*x4_4 = 1, b*x1_3*x4_5 = 1, z_aux*x3_0*c = 1, beta*x2_3 = 1, b*c*x1_0 = 1, -15*gama*x4_4*sgm = 1, -gama*sgm*x2_4*x4_0 = 1, x3_0*x4_1 = 1, x1_7*x4_1 = 1, b*x1_0*x4_5 = 1, -alpha*x1_0 = 1, x2_6 = 1, x2_7 = 1, x3_3*x4_2 = 1, delta*sgm*x3_5*x4_1 = 1, -x1_9 = 1, delta*sgm*x3_0*x4_6 = 1, b*x1_0*x4_4 = 1, x1_9*c = 1]
# 282, -5.577999884
# 2
# [x2_2 = 9, x3_2 = 14]
# [x2_2 = [4, 1, 4, 2, 2, 4, 4, 4, 4], x3_2 = [4, 2, 1, 4, 2, 2, 4, 2, 4, 2, 4, 2, 4, 2]]