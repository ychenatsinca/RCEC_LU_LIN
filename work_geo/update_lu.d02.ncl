begin 
basic = addfile("geo_em_2015.d02.nc","w")
RCEC = addfile("../d02.geo_RCEC_2015.nc","r")

;;;;URB_PARAM;;;;;;
;basic_URB_PARAM = basic->URB_PARAM
;RCEC_URB_PARAM = RCEC->URB_PARAM
;basic_URB_PARAM = (/RCEC_URB_PARAM/)
;basic->URB_PARAM = basic_URB_PARAM


;;;;LANDUSEF;;;;;
basic_LANDUSEF = basic->LANDUSEF
RCEC_LANDUSEF = RCEC->LANDUSEF
basic_LANDUSEF = (/RCEC_LANDUSEF/)
basic->LANDUSEF = basic_LANDUSEF

;;;;;LU_INDEX;;;;
basic_LU = basic->LU_INDEX
RCEC_LU = RCEC->LU_INDEX
basic_LU = (/RCEC_LU/)
basic->LU_INDEX = basic_LU

;;;;LANDMASK;;;;;;
;basic_LANDMASK = basic->LANDMASK
;RCEC_LANDMASK = RCEC->LANDMASK
;basic_LANDMASK = (/RCEC_LANDMASK/)
;basic->LANDMASK = basic_LANDMASK


end

