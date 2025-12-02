Sem - 1 project for university
Famous spinning torus with self implementation

This is a sample toy project to see if i could apply the theoretical math enough, here are some points i want to try later

1) increasing angular resolution (theta step and phi step to 0.03 and 0,01)
2) sampling multiple points around same pixel, similar to adding more augments, mapping each to brightness to reduce jaggedness you might see
3) Normalising lighting vector, non normalised ligting (rn) causes abrupt changes on every rotation, computinng L as dot product of unit normal and unit light diection changes this
4) smooth brightness mapping, interpolating (will have to learn) between two characters based on fractional index or just expanding brightness string to include gradations more


Open for input - 



edit 2: 

I wanted to try supersampling, but im not sure its practical given youd have to loop over pixel offsets , accumulate birhgtness per pixel, average and map to charset
