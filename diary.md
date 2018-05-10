## 日記

### 問題概要

2次元平面上に NC 個の町がある。
町と町の間に線分で道を敷設しすべての町を連結にしたい。
ただし道の敷設には長さに比例した費用がかかる。
道の端点は町か交差点である必要がある。
道の敷設の前に一度だけ、好きな数の交差点を建設してよい。
交差点の建設(の試行)には費用を支払う必要があり、ある確率で建設に失敗する。
できるだけ少ない費用を使ってすべての町を連結にせよ。

制約:

-   空間の幅 100 \le S \le 1000
-   町の数 10 \le NC \le 100
-   交差点の費用 0.0 \le junctionCost \le 10.0
-   交差点の建設の失敗の確率 0% \le failureProbability \le 40%

得点:

-   raw score: (道の総距離) + junctionCost * (交差点の数)
-   交差点の数には建設に失敗したものも含む
-   overall score: raw score が他の参加者より小さいと得点

### 5/10

今回は何すればいいか分からない。
焼き鈍しもビームサーチも効きにくそうでつらい。
chokudaiさんが1位取りそう。

下手に交差点の建設をさぼると連結化不可能になって詰むし難しいなあというのが第一印象だったが、visualizeしたら斜めの道もありだった。
道の交差の有無などは陽に確認したのに、道は軸平行だと思い込んでいた。
だって町の座標が整数点だったんだもん。
まだ数分考えただけとはいえ、誤読の情報だけ学習して他はちゃんと0から考察をrestartしないと前回と同じ轍を踏むぞ。

とりあえず交差点なしだと最小全域木やるだけ。
交差点ありだと(建設の成功失敗が決定した後でさえ)ちょっと面倒。
石鹸膜の表面積最小化みたいなあれが必要。
交差点ありの場合は最小Steiner木問題と呼ばれNP困難らしい (出典: [Steiner tree problem - Wikipedia](https://en.wikipedia.org/wiki/Steiner_tree_problem), [シュタイナー木 - Wikipedia](https://ja.wikipedia.org/wiki/%E3%82%B7%E3%83%A5%E3%82%BF%E3%82%A4%E3%83%8A%E3%83%BC%E6%9C%A8))。
でも今回は前回と違ってちゃんと平面グラフなので多項式解法あるかも。
そういえば師匠がICPC用としてSteiner木ライブラリ持ってたなあ。
交差点の位置が固定された後ですらNP困難なのに建設の試行からやるのつらいなあ。
Dreyfus-Wagner法というのがあるぽい (参考: [Spaghetti Source - 最小シュタイナー木 (Dreyfus-Wagner)](http://www.prefield.com/algorithm/dp/steiner_tree.html), [指数時間アルゴリズム入門, 岩田陽一 - SlideShare](https://www.slideshare.net/wata_orz/ss-12131479)) が、ほとんどの頂点が端点になる場合には実質ただの全探索になりそう。
Molle-Richter-Rossmanith法というのもあるらしいが方向性は変わらず。
すごく近いのがWikipediaにそのまま乗ってた: [Euclidean Steiner tree](https://en.wikipedia.org/wiki/Steiner_tree_problem#Euclidean_Steiner_tree)。

無駄なもの建設してもしかたがないので建設したものについては全部使うと仮定したい。
これは罠で、ふたつそろえば超効果的だけどひとつだけだとむしろ邪魔、みたいな場合が反例。
いやそんな例は実際に発生するのだろうか。
