infolevel[Groebner]:=10;
Et_hat := [52230949829341287542-w_0, c*q*w_0*y_0-c*w_0*x_0*y_0+b*w_0+w_1, 319346359494555263698-z_0, -c*q*w_0*y_0+h*z_0+z_1, 58508499888067911431-lm_0, lm_1, -w_1+19537856005373044245504078007367310677277386793046454350290085944259060292185008, ((-x_1*y_0+y_1*(q-x_0))*w_0+w_1*y_0*(q-x_0))*c+b*w_1+w_2, beta*v_0*x_0+d*x_0-lm_0+x_1, -beta*v_0*x_0+a*y_0+y_1, 45174929179497444858581413117111622188380755345680188803060076644450933095422168-z_1, -q*(w_0*y_1+w_1*y_0)*c+h*z_1+z_2, -w_2+7308460185655132347922186080241066725192061287343900955670854212756876199632368740786223137399622443121544920850579624473957733761708055000, ((q*y_2-x_0*y_2-2*x_1*y_1-x_2*y_0)*w_0+(-2*w_1*x_1+w_2*(q-x_0))*y_0+2*y_1*w_1*(q-x_0))*c+b*w_2+w_3, beta*x_0*v_1-lm_1+(beta*v_0+d)*x_1+x_2, -beta*v_0*x_1-beta*v_1*x_0+a*y_1+y_2, -k*y_0+u*v_0+v_1, 16898434055781315905253529446045954203513481637386942862788348479429461890719997505431007627452105248468135622571333248853228190140982965362-z_2, -q*(w_0*y_2+2*w_1*y_1+w_2*y_0)*c+h*z_2+z_3, -w_3+2733851159032861556295722638870108536845823965348298518825422166608634819109515940726390725102265226091835155656414924287017688228321788292565862217432956123400069400304908487612671854593565173153884, ((q*y_3-x_0*y_3-3*x_1*y_2-3*x_2*y_1-x_3*y_0)*w_0+(q*w_3-3*w_1*x_2-3*w_2*x_1-w_3*x_0)*y_0+(3*q*y_2-3*x_0*y_2-6*x_1*y_1)*w_1+3*y_1*w_2*(q-x_0))*c+b*w_3+w_4, (v_0*x_2+2*v_1*x_1+v_2*x_0)*beta+d*x_2-lm_2+x_3, (-v_0*x_2-2*v_1*x_1-v_2*x_0)*beta+a*y_2+y_3, lm_2, -k*y_1+u*v_1+v_2, 6321140480441290983083076174040379463683169888437081372250573798734496250742217770334479675984594512246544696026287548908924436324313336625570457859254090996399314409793513155252603063988505090799022-z_3, -3*(y_2*w_1+1/3*y_3*w_0+w_2*y_1+1/3*w_3*y_0)*q*c+h*z_3+z_4, -w_4+1022642522485788747927441717786096757195918261386291034871185077121310090986862490106281138432133559825983103366454074050515762595521877070112079319557185427424902363294627023129162791255866917512186147026097609576898165321452830864481050947445915935611806950, ((q*y_4-x_0*y_4-4*x_1*y_3-6*x_2*y_2-4*x_3*y_1-x_4*y_0)*w_0+(q*w_4-4*w_1*x_3-6*w_2*x_2-4*w_3*x_1-w_4*x_0)*y_0+(4*q*y_3-4*x_0*y_3-12*x_1*y_2-12*x_2*y_1)*w_1+(-6*w_2*y_2-4*w_3*y_1)*x_0+(4*q*w_3-12*w_2*x_1)*y_1+6*y_2*q*w_2)*c+b*w_4+w_5, (v_0*x_3+3*v_1*x_2+3*v_2*x_1+v_3*x_0)*beta+d*x_3+x_4-lm_3, (-v_0*x_3-3*v_1*x_2-3*v_2*x_1-v_3*x_0)*beta+a*y_3+y_4, lm_3, -k*y_2+u*v_2+v_3, 2364527792431954537422334788206631038607519238331768866218485936697035269328411184471851540505501258688559228294561234802309009195045303837743750753164170592800558427292618973209642264012440684462091989668715115185056359012930550670476152994430694392477025976-z_4, -4*q*(w_3*y_1+3/2*y_2*w_2+y_3*w_1+1/4*y_4*w_0+1/4*w_4*y_0)*c+h*z_4+z_5, -w_5+382536454239909189514623047661707636412443194727212981438558291786450379354533387234438820321006418542946991304925255215057428861681292651432751537557435728086670280072418971920890514934773221207369209850466172681630424720610682337190806662580608687082629545073330331069437890559155197190342836807496532332946322821568, ((q*y_5-x_0*y_5-5*x_1*y_4-10*x_2*y_3-10*x_3*y_2-5*x_4*y_1-x_5*y_0)*w_0+(q*w_5-5*w_1*x_4-10*w_2*x_3-10*w_3*x_2-5*w_4*x_1-w_5*x_0)*y_0+(5*q*y_4-5*x_0*y_4-20*x_1*y_3-30*x_2*y_2-20*x_3*y_1)*w_1+(-10*w_2*y_3-10*w_3*y_2-5*w_4*y_1)*x_0+(5*q*w_4-30*w_2*x_2-20*w_3*x_1)*y_1+(10*w_2*y_3+10*w_3*y_2)*q-30*x_1*y_2*w_2)*c+b*w_5+w_6, (v_0*x_4+4*v_1*x_3+6*v_2*x_2+4*v_3*x_1+v_4*x_0)*beta+d*x_4+x_5-lm_4, (-v_0*x_4-4*v_1*x_3-6*v_2*x_2-4*v_3*x_1-v_4*x_0)*beta+a*y_4+y_5, lm_4, -k*y_3+u*v_3+v_4, 884490970969975087779082219641565642058205060352732064564117765518601239958825298873030488541237926339345183294618351611610204686880251559653122416000044211775685227157647668787073817982674067517310718478363429427223184071477468799958456521796672262314261199751369386879783062566439036798134163182216100247040992527904-z_5, -q*(w_0*y_5+5*w_1*y_4+10*w_2*y_3+10*w_3*y_2+5*w_4*y_1+w_5*y_0)*c+h*z_5+z_6, -w_6+143094126837930000176400652442677244374272010662272372819977909007443081952821188439796836761061252408804932968701336217068786770127492197878765050026042371029260678587648743129845750982936711815327449874305375743728588861040951209044726644993997865326322711564394839093416748475135321489725237974647018404344025733213990051606094321485376718703702411361234489028805747750239104, ((q*y_6-x_0*y_6-6*x_1*y_5-15*x_2*y_4-20*x_3*y_3-15*x_4*y_2-6*x_5*y_1-x_6*y_0)*w_0+(q*w_6-6*w_1*x_5-15*w_2*x_4-20*w_3*x_3-15*w_4*x_2-6*w_5*x_1-w_6*x_0)*y_0+(6*q*y_5-6*x_0*y_5-30*x_1*y_4-60*x_2*y_3-60*x_3*y_2-30*x_4*y_1)*w_1+(-15*w_2*y_4-20*w_3*y_3-15*w_4*y_2-6*w_5*y_1)*x_0+(6*q*w_5-60*w_2*x_3-60*w_3*x_2-30*w_4*x_1)*y_1+(15*w_2*y_4+20*w_3*y_3+15*w_4*y_2)*q+(-60*x_1*y_3-90*x_2*y_2)*w_2-60*x_1*y_2*w_3)*c+b*w_6+w_7, (v_0*x_5+5*v_1*x_4+10*v_2*x_3+10*v_3*x_2+5*v_4*x_1+v_5*x_0)*beta+d*x_5+x_6-lm_5, (-v_0*x_5-5*v_1*x_4-10*v_2*x_3-10*v_3*x_2-5*v_4*x_1-v_5*x_0)*beta+a*y_5+y_6, lm_5, -k*y_4+u*v_4+v_5, -w_7+53526739500408612200960617219291534105334471311709617275186831957779114229820039339888962895202060020929635043761792015033747496835707384879577716086567363621559694383620886037780458960777921706810928418808898011425352262000503470998619485418201920345449771247349799660322549322867869865509168133628781606218246605075718634197432872344393452238681382505996025074407787779371130724351832627901171197630055931583152922552934146873727570804, 330858567292531723465007596962649273944618396796247859867224929090603383542643704601976910729119983490571888030154468843177597176996888753923920133400318235436721997202742461372251375895917081670151868026545977559254056851295584309402894139423004315280763626961374713520427194305665640131671672456227041842593096766788535094415958386535426346439253735248759858288385487808150498-z_6, -lm_1, -lm_2, -lm_3, -lm_4, -lm_5, -911586394381027768387213631862984339419230028117900591146545980600300772167418729010417021575335890456346812743330864860-v_3, z_aux-1];
vars:=[w_7, z_6, y_6, x_6, w_6, lm_5, z_5, y_5, x_5, w_5, v_5, lm_4, z_4, y_4, x_4, w_4, v_4, lm_3, z_3, y_3, x_3, w_3, v_3, lm_2, z_2, y_2, x_2, w_2, v_2, lm_1, z_1, y_1, x_1, w_1, v_1, lm_0, z_0, y_0, x_0, w_0, v_0, z_aux, w_aux, a, b, beta, c, d, h, k, q, u];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [lm, v_3] 22
# 257, [1 = 16, -1 = 2, v_4 = 1, a*y_0 = 1, b*w_1 = 1, -10*c*w_3*x_0*y_2 = 1, -w_0 = 1, -6*c*w_2*x_2*y_0 = 1, -2*c*w_1*x_0*y_1 = 1, -4*c*q*w_1*y_3 = 1, -15*c*w_4*x_2*y_0 = 1, -c*w_1*x_0*y_0 = 1, c*q*w_2*y_3 = 1, c*q*w_0*y_2 = 1, -5*c*w_1*x_0*y_4 = 1, -10*x_2*beta = 1, beta*v_0*x_5 = 1, -c*w_0*x_0*y_4 = 1, -beta*v_0*x_3 = 1, -c*w_0*x_5*y_0 = 1, c*q*w_1*y_2 = 1, u*v_4 = 1, -c*q*w_5*y_0 = 1, w_1 = 1, -k*y_3 = 1, -15*c*w_0*x_4*y_2 = 1, -3*c*w_1*x_2*y_0 = 1, -k*y_4 = 1, -c*w_0*x_1*y_0 = 1, -4*c*w_3*x_0*y_1 = 1, -2*c*w_1*x_1*y_0 = 1, -10*c*q*w_3*y_2 = 1, c*q*w_3*y_1 = 1, h*z_2 = 1, -6*c*w_0*x_5*y_1 = 1, -3*c*w_2*x_0*y_1 = 1, -6*c*w_0*x_2*y_2 = 1, b*w_6 = 1, -12*c*w_2*x_1*y_1 = 1, y_2 = 1, -c*w_0*x_2*y_0 = 1, -30*c*w_4*x_1*y_1 = 1, c*q*w_5*y_1 = 1, -c*w_3*x_0*y_0 = 1, -3*beta*v_1*x_2 = 1, c*q*w_2*y_4 = 1, -z_5 = 1, -15*c*w_4*x_0*y_2 = 1, -3*c*w_2*x_1*y_0 = 1, -beta*v_5*x_0 = 1, b*w_2 = 1, -c*w_0*x_0*y_1 = 1, v_1 = 1, -z_3 = 1, -w_6 = 1, -5*c*w_0*x_4*y_1 = 1, -c*q*w_0*y_4 = 1, beta*v_1*x_2 = 1, -5*c*q*w_4*y_1 = 1, -4*c*w_0*x_1*y_3 = 1, beta*v_0*x_1 = 1, -20*c*w_3*x_1*y_1 = 1, -c*w_0*x_0*y_0 = 1, -c*w_0*x_6*y_0 = 1, b*w_3 = 1, z_2 = 1, y_3 = 1, -z_0 = 1, c*q*w_5*y_0 = 1, x_6 = 1, -10*c*q*w_2*y_3 = 1, d*x_2 = 1, z_5 = 1, -c*q*w_0*y_2 = 1, -c*w_0*x_3*y_0 = 1, beta*v_2*x_2 = 1, -c*q*w_0*y_3 = 1, beta*v_1*x_1 = 1, -10*c*w_0*x_3*y_2 = 1, u*v_1 = 1, beta*v_1*x_4 = 1, -30*c*w_2*x_1*y_2 = 1, w_7 = 1, a*y_3 = 1, beta*x_0*v_1 = 1, -2*c*q*w_1*y_1 = 1, -w_2 = 1, -beta*x_0 = 1, z_4 = 1, -20*c*w_3*x_3*y_0 = 1, c*q*w_2*y_1 = 1, -4*beta*x_1 = 1, -c*q*w_4*y_0 = 1, beta*v_2*x_3 = 1, -6*c*w_1*x_0*y_5 = 1, -w_5 = 1, -6*c*w_2*x_0*y_2 = 1, y_1 = 1, -6*c*q*w_2*y_2 = 1, c*q*w_2*y_2 = 1, -w_7 = 1, -10*c*w_0*x_2*y_3 = 1, -3*c*w_1*x_0*y_2 = 1, -10*c*w_2*x_3*y_0 = 1, -12*c*w_1*x_2*y_1 = 1, d*x_4 = 1, beta*v_2*x_1 = 1, -60*c*w_3*x_2*y_1 = 1, c*q*w_6*y_0 = 1, v_2 = 1, -c*q*w_0*y_0 = 1, c*q*w_1*y_5 = 1, -c*q*w_0*y_5 = 1, -c*q*w_1*y_0 = 1, a*y_1 = 1, -4*c*w_1*x_3*y_0 = 1, beta*v_2*x_0 = 1, -c*q*w_0*y_1 = 1, c*q*w_0*y_1 = 1, -3*c*q*w_1*y_2 = 1, beta*v_0*x_2 = 1, -30*c*w_1*x_1*y_4 = 1, x_3 = 1, -5*c*w_0*x_1*y_4 = 1, -w_3 = 1, -5*c*q*w_1*y_4 = 1, -c*w_2*x_0*y_0 = 1, -4*c*q*w_3*y_1 = 1, -30*c*w_2*x_2*y_1 = 1, h*z_4 = 1, -c*q*w_2*y_0 = 1, -6*c*w_1*x_1*y_1 = 1, -20*c*w_1*x_3*y_1 = 1, -3*c*q*w_2*y_1 = 1, -beta*v_0*x_4 = 1, c*q*w_0*y_4 = 1, -3*c*w_0*x_1*y_2 = 1, -20*c*w_3*x_0*y_3 = 1, -5*beta*v_1*x_4 = 1, b*w_4 = 1, -z_6 = 1, c*q*w_4*y_2 = 1, -4*c*w_1*x_0*y_3 = 1, v_5 = 1, h*z_0 = 1, beta*v_0*x_3 = 1, -beta*v_0*x_1 = 1, c*q*w_0*y_6 = 1, -c*w_0*x_0*y_5 = 1, -beta*v_0*x_2 = 1, -c*w_0*x_0*y_6 = 1, c*q*w_4*y_0 = 1, x_2*beta = 1, -60*c*w_1*x_2*y_3 = 1, c*q*w_3*y_0 = 1, -12*c*w_1*x_1*y_2 = 1, beta*v_5*x_0 = 1, c*q*w_0*y_3 = 1, b*w_5 = 1, -c*w_0*x_0*y_3 = 1, x_2 = 1, beta*v_1*x_3 = 1, -5*c*w_1*x_4*y_0 = 1, -c*q*w_3*y_0 = 1, z_1 = 1, -5*c*w_4*x_1*y_0 = 1, -60*c*w_3*x_1*y_2 = 1, beta*v_0*x_4 = 1, -w_4 = 1, -w_1 = 1, c*q*w_3*y_2 = 1, -beta*x_0*v_1 = 1, beta*v_0*x_0 = 1, -30*c*w_1*x_4*y_1 = 1, x_1 = 1, beta*v_4*x_0 = 1, -10*c*w_3*x_2*y_0 = 1, -15*c*w_2*x_4*y_0 = 1, c*q*w_1*y_1 = 1, c*q*w_0*y_5 = 1, -4*beta*v_1*x_3 = 1, -3*c*w_0*x_2*y_1 = 1, -5*c*w_4*x_0*y_1 = 1, -z_1 = 1, beta*x_1 = 1, -z_4 = 1, y_6 = 1, h*z_5 = 1, -10*c*w_2*x_0*y_3 = 1, -6*beta*v_2*x_2 = 1, c*q*w_0*y_0 = 1, x_4 = 1, -beta*v_0*x_5 = 1, -60*c*w_1*x_3*y_2 = 1, w_5 = 1, w_2 = 1, d*x_1 = 1, -6*c*w_5*x_0*y_1 = 1, -6*c*w_1*x_5*y_0 = 1, c*q*w_1*y_3 = 1, h*z_1 = 1, h*z_3 = 1, w_3 = 1, c*q*w_3*y_3 = 1, d*x_5 = 1, -beta*v_0*x_0 = 1, -2*c*w_0*x_1*y_1 = 1, -20*c*w_0*x_3*y_3 = 1, x_5 = 1, b*w_0 = 1, -k*y_0 = 1, -4*c*w_0*x_3*y_1 = 1, -15*c*w_0*x_2*y_4 = 1, beta*v_4*x_1 = 1, -c*w_4*x_0*y_0 = 1, -5*beta*v_4*x_1 = 1, y_5 = 1, a*y_2 = 1, -3*beta*v_2*x_1 = 1, u = 1, -20*c*w_1*x_1*y_3 = 1, y_4 = 1, w_6 = 1, -60*c*w_2*x_3*y_1 = 1, -k*y_1 = 1, u*v_0 = 1, -10*beta*v_2*x_3 = 1, u*v_2 = 1, -c*w_0*x_0*y_2 = 1, a*y_5 = 1, -30*c*w_1*x_2*y_2 = 1, w_4 = 1, a*y_4 = 1, c*q*w_4*y_1 = 1, -z_2 = 1, -60*c*w_2*x_1*y_3 = 1, -2*beta*v_1*x_1 = 1, -c*w_5*x_0*y_0 = 1, -4*c*w_3*x_1*y_0 = 1, -6*c*w_5*x_1*y_0 = 1, -beta*v_4*x_0 = 1, -c*w_0*x_4*y_0 = 1, -k*y_2 = 1, d*x_0 = 1, z_3 = 1, d*x_3 = 1, -90*c*w_2*x_2*y_2 = 1, z_6 = 1, -beta*v_2*x_0 = 1, c*q*w_1*y_0 = 1, -c*w_6*x_0*y_0 = 1, c*q*w_1*y_4 = 1, -6*c*w_0*x_1*y_5 = 1, c*q*w_2*y_0 = 1, beta*x_0 = 1, z_aux = 1, -15*c*w_2*x_0*y_4 = 1]
# 273, -5.441897751
# 2
# [v_3 = 8, lm = 1]
# [v_3 = [3, 3, 1, 3, 3, 2, 3, 3], lm = [1]]
# [v_3 = [1, 1, 1, 1, 1, 1, 1, 1], lm = [1]]