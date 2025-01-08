#pragma once
#include <random>
//#include <generator>

class discreteDistributions {
public:
  enum Distribution{ uniform, mito, numt };
  std::discrete_distribution<> d;
  discreteDistributions(Distribution dist_p = uniform) {

    switch(dist_p){
    case mito:
      d = {
		12221, 12807, 13424, 14074, 14759, 15479, 16239, 17038, 17881, 18768, 
19702, 20687, 21724, 22816, 23966, 25178, 26455, 27799, 29215, 30706, 
32276, 33929, 35670, 37502, 39431, 41461, 43597, 45845, 48210, 50696, 
53312, 56061, 58951, 61988, 65179, 68530, 72049, 75743, 79619, 83686, 
 87950,  92420,  97104, 102009, 107146, 112521, 118144, 124022, 130165,136581, 
143278, 150265, 157549, 165140, 173045, 181270, 189825, 198715, 207946, 217526, 
 227459, 237751, 248405, 259426, 270817, 282579, 294714, 307223, 320105, 333360,
 346985, 360977, 375332, 390045, 405110, 420521, 436269, 452346, 468742,485447, 
 502448, 519734, 537293, 555109, 573170, 591459, 609962, 628662, 647544, 666591, 
 685785, 705110, 724547, 744080, 763691, 783362, 803075, 822813, 842559, 862294,
   882001,  901663,  921263,  940785,  960211,  979524,  998709, 1017750, 1036630, 1055333, 
  1073845, 1092149, 1110231, 1128075, 1145668, 1162994, 1180039, 1196790, 1213233, 1229354, 
  1245141, 1260581, 1275662, 1290372, 1304699, 1318635, 1332167, 1345288, 1357987, 1370257, 
  1382090, 1393478, 1404416, 1414898, 1424918, 1434470, 1443551, 1452155, 1460278, 1467916, 
  1475064, 1481718, 1487872, 1493520, 1498658, 1503278, 1507373, 1510936, 1513959, 1516433, 
  1518349, 1519700, 1520476, 1520669, 1520272, 1519277, 1517681, 1515480, 1512673, 1509260, 
  1505245, 1500636, 1495441, 1489674, 1483350, 1476489, 1469111, 1461241, 1452905, 1444130, 
  1434945, 1425381, 1415467, 1405233, 1394709, 1383924, 1372904, 1361677, 1350267, 1338698, 
  1326990, 1315165, 1303239, 1291232, 1279156, 1267027, 1254856, 1242656, 1230436, 1218206, 
  1205974, 1193747, 1181532, 1169336, 1157164, 1145021, 1132912, 1120841, 1108812, 1096829, 
  1084895, 1073013, 1061187, 1049419, 1037711, 1026066, 1014487, 1002976,991534,  980164, 
  968867, 957645, 946501, 935435, 924448, 913544, 902722, 891983, 881330,870764, 
  860284, 849893, 839591, 829379, 819257, 809227, 799289, 789444, 779691,770033, 
  760468, 750997, 741622, 732341, 723155, 714064, 705068, 696167, 687362, 678652, 
  670037, 661517, 653092, 644761, 636524, 628382, 620333, 612377, 604514,596744,
  589066, 581479, 573984, 566578, 559263, 552038, 544901, 537852, 530891, 524017, 
  517230, 510528, 503911, 497378, 490928, 484562, 478277, 472074, 465951, 459908, 
  453944, 448058, 442250, 436519, 430863, 425283, 419776, 414344, 408984, 403696, 
  398479, 393332, 388255, 383246, 378306, 373432, 368625, 363883, 359206, 354593, 
  350043, 345556, 341130, 336765, 332460, 328214, 324027, 319897, 315825, 311809,
   307848, 303942, 300090, 296292, 292546, 288853, 285210, 281618, 278076, 274584, 
   271139, 267743, 264394, 261091, 257834, 254623, 251456, 248333, 245253,
 242217, 239222, 236269, 233357, 230486, 227654, 224861, 222108, 219392, 216714,
 214073, 211469, 208900, 206368, 203870, 201406, 198977, 196581, 194218, 191887,
 189589, 187322, 185086, 182881, 180706, 178561, 176446, 174359, 172301, 170270,
 168268, 166292, 164344, 162422, 160526, 158656, 156811, 154991, 153196, 151424,
 149677, 147953, 146253, 144575, 142920, 141287, 139676, 138086, 136517, 134970,
 133443, 131936, 130449, 128982, 127535, 126106, 124697, 123306, 121933, 120579,
 119242, 117922, 116620, 115335, 114067, 112815, 111580, 110360, 109157, 107969,
 106796, 105638, 104496, 103368, 102254, 101155, 100070,  98999,  97942,  96898,
  95867, 94849, 93845, 92853, 91873, 90906, 89951, 89009, 88078, 87158, 86250, 
  85354, 84469, 83594, 82731, 81878, 81036, 80205, 79383, 78572, 77771, 
  76979, 76198, 75426, 74663, 73910, 73165, 72430, 71704, 70987, 70278, 
  69578, 68886, 68203, 67528, 66861, 66201, 65550, 64907, 64271, 63643, 
  63022, 62409, 61803, 61204, 60612, 60027, 59449, 58877, 58313, 57755, 
  57203, 56658, 56119, 55587, 55060, 54540, 54026, 53517, 53015, 52518, 
  52027, 51542, 51062, 50587, 50118, 49654, 49196, 48743, 48294, 47851, 
  47413, 46980, 46551, 46128, 45709, 45295, 44885, 44480, 44080, 43684, 
  43292, 42904, 42521, 42142, 41767, 41397, 41030, 40667, 40309, 39954, 
  39603, 39256, 38912, 38573, 38237, 37904, 37575, 37250, 36928, 36609,
};
      break;

      
    case numt :
      d = {
14, 15, 15, 16, 17, 18, 19, 20, 21, 22, 
23, 24, 25, 26, 28, 29, 31, 32, 34, 36, 
38, 40, 42, 44, 47, 49, 52, 55, 58, 61, 
 65,  69,  73,  77,  81,  86,  91,  96, 102, 108,
 114, 121, 128, 136, 144, 153, 162, 172, 183, 194, 
 206, 219, 233, 248, 264, 281, 299, 318, 338, 360,
  384, 409, 436, 465, 496, 529, 565, 603, 644, 688, 
   735,  786,  841,  899,  962, 1030, 1103, 1181, 1265, 1356, 
  1454, 1560, 1674, 1796, 1929, 2072, 2227, 2394, 2574, 2769, 
  2981, 3209, 3456, 3723, 4013, 4326, 4665, 5033, 5431, 5862,
    6328,  6833,  7380,  7972,  8613,  9306, 10056, 10867, 11744, 12691, 
   13714, 14819, 16011, 17297, 18684, 20178, 21787, 23519, 25383, 27387, 
   29541, 31854, 34338, 37003, 39861, 42925, 46207, 49723, 53487, 57515, 
    61824,  66432,  71359,  76625,  82252,  88263,  94685, 101544, 108868, 116690, 
   125043, 133961, 143483, 153651, 164509, 176102, 188483, 201703, 215821, 230899, 
   247000, 264195, 282557, 302162, 323093, 345434, 369274, 394702, 421808, 450682, 
   481413, 514047, 548667, 584659, 619537, 649576, 664393, 657050, 634984, 603363,
    566222, 530450, 496431, 464234, 433847, 405227, 378310, 353026, 329301, 307058, 
 286219, 266711, 248458, 231390, 215439, 200539, 186628, 173645, 161534, 150240, 
 139714, 129906, 120770, 112265, 104348,  96981,  90130,  83758,  77836, 72332,
  67218, 62469, 58059, 53965, 50166, 46642, 43372, 40341, 37530, 34926, 
  32512, 30277, 28207, 26291, 24518, 22879, 21362, 19961, 18667, 17472, 
  16370, 15354, 14417, 13556, 12763, 12035, 11367, 10756, 10196,  9685, 
  9220, 8797, 8414, 8067, 7755, 7476, 7227, 7006, 6812, 6643, 
  6498, 6375, 6273, 6190, 6127, 6080, 6051, 6037, 6039, 6055, 
  6084, 6127, 6182, 6250, 6329, 6419, 6521, 6633, 6756, 6888, 
  7031, 7183, 7344, 7515, 7695, 7884, 8082, 8288, 8504, 8727, 
   8959,  9200,  9449,  9706,  9971, 10244, 10526, 10815, 11112, 11417, 
  11730, 12050, 12378, 12714, 13057, 13407, 13764, 14128, 14500, 14877, 
  15262, 15652, 16049, 16452, 16860, 17273, 17691, 18114, 18540, 18969, 
  19399, 19826, 20252, 20675, 21090, 21491, 21873, 22195, 22480, 22746, 
  23001, 23230, 23423, 23582, 23683, 23740, 23765, 23763, 23730, 23663, 
  23567, 23439, 23268, 23069, 22851, 22619, 22364, 22090, 21800, 21495, 
  21181, 20859, 20529, 20184, 19829, 19463, 19092, 18716, 18339, 17961, 
  17583, 17205, 16829, 16453, 16080, 15708, 15339, 14973, 14611, 14251, 
  13894, 13541, 13193, 12848, 12508, 12172, 11841, 11515, 11194, 10877, 
  10566, 10261,  9960,  9665,  9376,  9092,  8814,  8542,  8275,  8014, 
  7759, 7510, 7266, 7029, 6797, 6571, 6351, 6137, 5929, 5726, 
  5529, 5337, 5151, 4970, 4794, 4624, 4458, 4298, 4142, 3992, 
  3846, 3704, 3567, 3435, 3306, 3182, 3062, 2946, 2834, 2726, 
  2621, 2520, 2422, 2328, 2237, 2150, 2065, 1984, 1905, 1829, 
  1756, 1686, 1618, 1553, 1491, 1430, 1372, 1316, 1262, 1211, 
  1161, 1113, 1067, 1023,  980,  940,  900,  863,  827,  792, 
  759, 727, 696, 667, 639, 612, 586, 561, 537, 514, 
  492, 471, 451, 432, 413, 395, 378, 362, 347, 332, 
  317, 304, 290, 278, 266, 254, 243, 233, 223, 213,
  204, 195, 186, 178, 170, 163, 156, 149, 142, 136, 
  130, 124, 119, 114, 109, 104,  99,  95,  91,  87, 
  83, 79, 76, 72, 69, 66, 63, 60, 58, 55, 
  53, 50, 48, 46, 44, 42, 40, 38, 37, 35 
};
      break;
    case uniform:
       d = {
      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
};
  break;

  }
  }

  template<class Generator>
  int draw(Generator& gen) {
    return d(gen) + 1;
      }
};