infolevel[Groebner]:=10;
Et_hat := [67651370044379-Tr_0, -In_0*g+Tr_0*nu+Tr_1, 44197926932887-b_0, b_1, -Tr_1+595081832352031773732919959, -In_1*g+Tr_1*nu+Tr_2, -S_0*Tr_0*b_0*d+In_0*N*a+In_0*N*g-In_0*S_0*b_0+In_1*N, -Tr_2+290576211942135810635084372010790423432995172706338758857160810675143/10292248380901, -In_2*g+Tr_2*nu+Tr_3, ((-Tr_1*d-In_1)*b_0-b_1*(Tr_0*d+In_0))*S_0-S_1*(Tr_0*d+In_0)*b_0+N*((a+g)*In_1+In_2), S_0*Tr_0*b_0*d+In_0*S_0*b_0+N*S_1, -Tr_3-70745425813472644718004806565247683207475770633280399186993167536549312315451607218330468368290394642040469283/317791130202477767944715403, -In_3*g+Tr_3*nu+Tr_4, ((-Tr_2*d-In_2)*b_0+(-Tr_0*b_2-2*Tr_1*b_1)*d-In_0*b_2-2*In_1*b_1)*S_0+((-2*S_1*Tr_1-S_2*Tr_0)*d-S_2*In_0-2*S_1*In_1)*b_0-2*S_1*Tr_0*b_1*d-2*In_0*S_1*b_1+((a+g)*In_2+In_3)*N, ((Tr_1*d+In_1)*b_0+b_1*(Tr_0*d+In_0))*S_0+S_1*(Tr_0*d+In_0)*b_0+N*S_2, b_2, -Tr_4+17224105304690271905445189600002519897919397569284069169650678285588362178462316101542095209616100137137742844747871808922468025374235265181172324855431/9812355735873452062282336515017957154309, -In_4*g+Tr_4*nu+Tr_5, ((-Tr_0*b_3-3*Tr_1*b_2-3*Tr_2*b_1-Tr_3*b_0)*S_0+(-3*S_1*Tr_2-3*S_2*Tr_1-S_3*Tr_0)*b_0+(-3*Tr_0*b_2-6*Tr_1*b_1)*S_1-3*Tr_0*S_2*b_1)*d+(-In_0*b_3-3*In_1*b_2-3*In_2*b_1-In_3*b_0)*S_0+(-In_0*S_3-3*In_1*S_2-3*In_2*S_1)*b_0+(-3*In_0*b_2-6*In_1*b_1)*S_1-3*In_0*S_2*b_1+((a+g)*In_3+In_4)*N, ((Tr_2*d+In_2)*b_0+(Tr_0*b_2+2*Tr_1*b_1)*d+In_0*b_2+2*In_1*b_1)*S_0+((2*S_1*Tr_1+S_2*Tr_0)*d+S_2*In_0+2*S_1*In_1)*b_0+2*S_1*Tr_0*b_1*d+2*In_0*S_1*b_1+N*S_3, b_3, -Tr_5-4193483891513319245117748297792667189570199919863868695152186462738827486131555191870270443883373541735808921585902267948221907336992032959958350207437105902464854474367316529587851402962960771/302973607306104532173165052222274401177132183996357227, -In_5*g+Tr_5*nu+Tr_6, ((-Tr_0*b_4-4*Tr_1*b_3-6*Tr_2*b_2-4*Tr_3*b_1-Tr_4*b_0)*S_0+(-4*S_1*Tr_3-6*S_2*Tr_2-4*S_3*Tr_1-S_4*Tr_0)*b_0+(-4*Tr_0*b_3-12*Tr_1*b_2-12*Tr_2*b_1)*S_1+(-12*S_2*Tr_1-4*S_3*Tr_0)*b_1-6*b_2*Tr_0*S_2)*d+(-In_0*b_4-4*In_1*b_3-6*In_2*b_2-4*In_3*b_1-In_4*b_0)*S_0+(-In_0*S_4-4*In_1*S_3-6*In_2*S_2-4*In_3*S_1)*b_0+(-4*In_0*b_3-12*In_1*b_2-12*In_2*b_1)*S_1+(-4*In_0*S_3-12*In_1*S_2)*b_1-6*In_0*S_2*b_2+N*((a+g)*In_4+In_5), ((Tr_0*b_3+3*Tr_1*b_2+3*Tr_2*b_1+Tr_3*b_0)*S_0+(3*S_1*Tr_2+3*S_2*Tr_1+S_3*Tr_0)*b_0+(3*Tr_0*b_2+6*Tr_1*b_1)*S_1+3*Tr_0*S_2*b_1)*d+(In_0*b_3+3*In_1*b_2+3*In_2*b_1+In_3*b_0)*S_0+(In_0*S_3+3*In_1*S_2+3*In_2*S_1)*b_0+(3*In_0*b_2+6*In_1*b_1)*S_1+3*In_0*S_2*b_1+N*S_4, b_4, -Tr_6+1020970717334841161026087834860634151526464762196980998307424930001934521573700943261836088852847824851132114739414525895015805914753416210481659580547218425355790507901202624710432379958121013760669637201485567392234084105307644924327/9354838857755969266658148215086022452449234147728897079818280364581, -In_6*g+Tr_6*nu+Tr_7, ((-Tr_0*b_5-5*Tr_1*b_4-10*Tr_2*b_3-10*Tr_3*b_2-5*Tr_4*b_1-Tr_5*b_0)*S_0+(-5*S_1*Tr_4-10*S_2*Tr_3-10*S_3*Tr_2-5*S_4*Tr_1-S_5*Tr_0)*b_0+(-5*Tr_0*b_4-20*Tr_1*b_3-30*Tr_2*b_2-20*Tr_3*b_1)*S_1+(-30*S_2*Tr_2-20*S_3*Tr_1-5*S_4*Tr_0)*b_1+(-10*Tr_0*b_3-30*Tr_1*b_2)*S_2-10*b_2*Tr_0*S_3)*d+(-In_0*b_5-5*In_1*b_4-10*In_2*b_3-10*In_3*b_2-5*In_4*b_1-In_5*b_0)*S_0+(-In_0*S_5-5*In_1*S_4-10*In_2*S_3-10*In_3*S_2-5*In_4*S_1)*b_0+(-5*In_0*b_4-20*In_1*b_3-30*In_2*b_2-20*In_3*b_1)*S_1+(-5*In_0*S_4-20*In_1*S_3-30*In_2*S_2)*b_1+(-10*In_0*b_3-30*In_1*b_2)*S_2-10*In_0*S_3*b_2+((a+g)*In_5+In_6)*N, ((Tr_0*b_4+4*Tr_1*b_3+6*Tr_2*b_2+4*Tr_3*b_1+Tr_4*b_0)*S_0+(4*S_1*Tr_3+6*S_2*Tr_2+4*S_3*Tr_1+S_4*Tr_0)*b_0+(4*Tr_0*b_3+12*Tr_1*b_2+12*Tr_2*b_1)*S_1+(12*S_2*Tr_1+4*S_3*Tr_0)*b_1+6*b_2*Tr_0*S_2)*d+(In_0*b_4+4*In_1*b_3+6*In_2*b_2+4*In_3*b_1+In_4*b_0)*S_0+(In_0*S_4+4*In_1*S_3+6*In_2*S_2+4*In_3*S_1)*b_0+(4*In_0*b_3+12*In_1*b_2+12*In_2*b_1)*S_1+(4*In_0*S_3+12*In_1*S_2)*b_1+6*In_0*S_2*b_2+N*S_5, b_5, -Tr_7-248571648925284210780818375610177694319782421054692333704901283512284387618156589053958573434984085246793236214742575370030818724761592268461121242318681902847893159611021059907262039522907561783255817856170152756415637278879515767574516916721723252798274485193956438270237027/288846975261985904792790726869333993065970290256580544848150489136982631711802443, -b_1, -b_2, -b_3, -b_4, -b_5, 1569521626332760647736640829391165386490509916207254062270790234169962922631004625529238478747280/317791130202477767944715403-S_2, N*z_aux-1];
vars:=[Tr_7, Tr_6, In_6, Tr_5, In_5, b_5, S_5, Tr_4, In_4, b_4, S_4, Tr_3, In_3, b_3, S_3, Tr_2, In_2, b_2, S_2, Tr_1, In_1, b_1, S_1, Tr_0, In_0, b_0, S_0, z_aux, w_aux, N, a, d, g, nu];
gb:=CodeTools[Usage](Groebner[Basis](Et_hat, tdeg(op(vars)), characteristic=11863279), output='all');
# [S_2, b] 254
# 131, [1 = 4, -Tr_4 = 1, -10*In_3 = 1, -6*In_2 = 1, S_1*In_1 = 1, -3*S_1*Tr_2*d = 1, S_4*Tr_0*d = 1, In_2 = 1, -179309634593493733658162212844461858518738205059942771396968370122660798111006850172798173061927452435076970925532086626741233547613316281809745770255690419192793534461037590319412407502817305399172689471404860699487289686616281977054599048429017935279549578831136/246035506446721274012617103349915219301813950982064232570241051996976213527 = 1, In_0 = 1, -S_0*Tr_5*d = 1, -195729598087814495425614884200452077342696811229758649392798447939034136730600123569559915378863964696973980549415713899267784372098356115932139834648691705675571488195208248627151963570/274850520894250566433945179777030662974823983471063 = 1, N*In_2 = 1, -5*S_4*Tr_1*d = 1, -In_2*S_0 = 1, -Tr_3 = 1, N*S_1 = 1, In_0*S_0 = 1, -In_4*S_0 = 1, Tr_1*d = 1, Tr_5 = 1, -10*S_3*Tr_2*d = 1, Tr_2*d = 1, -Tr_6 = 1, S_3*Tr_0*d = 1, S_1*Tr_2*d = 1, In_3*S_0 = 1, In_2*S_0 = 1, S_0*Tr_4*d = 1, In_1*N*a = 1, -10*Tr_3*d = 1, -Tr_7 = 1, -g*In_1 = 1, Tr_7 = 1, -In_0 = 1, -In_0*S_1 = 1, N*S_5 = 1, In_0*S_4 = 1, Tr_3 = 1, -In_3*S_0 = 1, S_1*Tr_3*d = 1, N*In_4 = 1, Tr_4 = 1, In_4*N*a = 1, In_1 = 1, -3*In_2*S_1 = 1, In_2*N*g = 1, -5*S_1*Tr_4*d = 1, nu*Tr_5 = 1, -Tr_1*d*S_0 = 1, Tr_0*nu = 1, nu*Tr_3 = 1, -4*In_3*S_1 = 1, -S_3*In_0 = 1, -In_5*g = 1, nu*Tr_4 = 1, -4360269199474597815095572546544347080173935020453374789825453251201657350086524463723247049405869120527092/6266128001226122488850503 = 1, -In_5*S_0 = 1, Tr_6 = 1, -S_0*Tr_2*d = 1, -6*Tr_2*d = 1, S_3*In_0 = 1, In_1*S_3 = 1, nu*Tr_2 = 1, -g*In_4 = 1, -Tr_2 = 1, In_1*N*g = 1, Tr_1*d*S_0 = 1, -Tr_0 = 1, In_5*N*a = 1, -4*S_1*Tr_3*d = 1, S_0*Tr_2*d = 1, Tr_0*S_0*d = 1, -S_1*Tr_0*d = 1, S_0*In_1 = 1, In_4*N*g = 1, S_1*Tr_0*d = 1, In_2*N*a = 1, -S_0*Tr_3*d = 1, -3*In_1 = 1, Tr_2 = 1, In_3*N*g = 1, N*In_6 = 1, N = 1, nu*Tr_6 = 1, N*In_5 = 1, S_1*Tr_1*d = 1, -4*S_3*Tr_1*d = 1, -Tr_0*S_0*d = 1, In_0*N*g = 1, -3*Tr_1*d = 1, N*In_1 = 1, -2*S_1*In_1 = 1, -10*In_2*S_3 = 1, -g*In_2 = 1, In_0*N*a = 1, -In_0*g = 1, N*S_3 = 1, z_aux*N = 1, N*S_4 = 1, -S_0*Tr_4*d = 1, Tr_1 = 1, In_3*N*a = 1, S_0*Tr_3*d = 1, -g*In_3 = 1, In_5*N*g = 1, -S_3*Tr_0*d = 1, -5429305681874946826554610 = 1, -Tr_0*d = 1, In_2*S_1 = 1, -2*S_1*Tr_1*d = 1, nu*Tr_1 = 1, -S_0*In_1 = 1, S_3*Tr_1*d = 1, -4*In_1*S_3 = 1, -5*In_1*S_4 = 1, In_4*S_0 = 1, -5*In_4*S_1 = 1, -In_0*S_5 = 1, N*In_3 = 1, -In_0*S_4 = 1, In_0*S_1 = 1, -In_6*g = 1, In_3*S_1 = 1, -Tr_5 = 1, -Tr_1 = 1, Tr_0*d = 1, -S_4*Tr_0*d = 1, -S_5*Tr_0*d = 1, -In_0*S_0 = 1, -1 = 1]
# 134, -4.856457879
# 1
# [S_2 = 15, b = 58]
# [S_2 = [4, 3, 2, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3, 4, 3], b = [4, 3, 4, 4, 3, 3, 4, 3, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3]]