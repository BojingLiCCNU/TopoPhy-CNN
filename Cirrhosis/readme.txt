cirr_raw_all.txt
	第一列是二进制0或1的标签，第二列是乘的倍数（用于NB分布合成数据），第三列为标签，之后数据为丰度
cirr_raw_data.txt
	将PopPhy-CNN中丰度数据转为NB分布可以合成的数据，过程如下：1，删除微生物特征的名字即第一列
	2，转置行列，使得每一行为一个样本 3，数据同乘100000，保证NB的输入数据为整数
cirr_raw_label.txt
	将PopPhy-CNN的label数据转换为0或1的二进制标签数据，作为NB分布的输入


nb_cirr_trainabun_10.txt
	使用NB分布，对cirr数据，扩充为原来的10倍后的丰度数据，可直接用于PopPhy-CNN训练
	如何得到：augment/results文件夹中，找到nb_cirr_data_10.txt文件，改文件为cirr数据集
	经过NB分布扩充10倍后的数据，将nb_cirr_data_10.txt文件中的数据转置，除100000，并加上微生物
	特征的名字即第一列

nb_cirr_trainlabel_10.txt
	使用NB分布，对cirr数据，扩充为原来的10倍后的label数据，可直接用于PopPhy-CNN训练
	如何得到：augment/results文件夹中，找到nb_cirr_label_10.txt文件，将0或1标签改为n或cirrhosis标签