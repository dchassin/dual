[![Test](https://github.com/dchassin/dual/actions/workflows/test.yml/badge.svg)](https://github.com/dchassin/dual/actions/workflows/test.yml)

# Dual number arithmetic

For information on dual numbers see https://www.researchgate.net/profile/Farid-Messelmi/publication/264300733_Analysis_of_Dual_Functions/links/53d7c11e0cf2a19eee7fcf1c/Analysis-of-Dual-Functions.pdf.

## Installation

~~~
python3 -m venv .
. bin/activate
python3 -m pip install pip --upgrade -r requirements.txt
~~~

## Example

The following commands

~~~
form dual import dual
x = dual(1,2)
y = dual(3,4)
x*y
~~~

generates the following output

~~~
<dual:3+10w>
~~~
