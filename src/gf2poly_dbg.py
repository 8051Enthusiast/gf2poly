import gdb

obj = gdb.current_objfile()
if obj:
    limbs_sym = obj.lookup_static_symbol("gf2poly::Gf2Poly::limbs")
else:
    limbs_sym = None


class Gf2PolyPrinter:
    def __init__(self, val):
        self.val = val
        if limbs_sym:
            self.limbs_fun = limbs_sym.value()
    
    def to_string_brute(self):
        try:
            limbs = self.val["limbs"]
            word_type = limbs.type.template_argument(0)
            word_bits = word_type.sizeof * 8
            arr_len = limbs["len"]
            arr_ptr = limbs["buf"]["inner"]["ptr"]["pointer"]["pointer"].cast(word_type.pointer())
            arr = [arr_ptr[i] for i in range(int(arr_len))]
            return to_string_from_list(arr, word_bits)
        except gdb.error:
            return None
    
    def to_string_sym(self):
        try:
            arr = self.limbs_fun(self.val.address).to_array()
            if not arr:
                return None
            size = arr.type.sizeof
            if size == 0:
                word_width = 1
            else:
                word_width = arr[0].type.sizeof
            arr = [arr[i] for i in range(size // word_width)]
            return to_string_from_list(arr, word_width * 8)
        except gdb.error:
            return None
        
    def to_string(self):
        ret = None
        ret = self.to_string_brute()
        if not ret and self.limbs_fun:
            ret = self.to_string_sym()
        return ret
 


def to_string_from_list(arr, bits_per_word):
    if len(arr) == 0:
        return "Gf2Poly(0)"
    digits = bits_per_word // 4
    front = arr[-1]
    ret = f"{int(front):x}"
    for word in arr[-2::-1]:
        ret += f"{int(word):0{digits}x}"
    return f"Gf2Poly({ret})"


def lookup(val):
    lookup_tag = val.type.tag
    if lookup_tag == None:
        return None
    if lookup_tag == "gf2poly::Gf2Poly":
        return Gf2PolyPrinter(val)
    return None

if obj:
    obj.pretty_printers.append(lookup)