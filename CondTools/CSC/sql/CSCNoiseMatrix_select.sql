REM For CSCNOISEMATRIX:
SELECT
 record_id iov_value_id,
 run_num time
FROM NOISEMATRIX;

REM For CSCNOISEMATRIX_MAP:
SELECT
 map_index map_id,
 record_id iov_value_id,
 layer_id csc_int_id
FROM NOISEMATRIX_MAP
 order by iov_value_id,map_id;

REM For CSCNOISEMATRIX_DATA:
SELECT
 NOISEMATRIX_DATA.vec_index vec_index,
 NOISEMATRIX_MAP.map_index map_id,
 NOISEMATRIX_MAP.record_id iov_value_id,
 NOISEMATRIX_DATA.elem33 cscmatrix_elem33,
 NOISEMATRIX_DATA.elem34 cscmatrix_elem34,
 NOISEMATRIX_DATA.elem35 cscmatrix_elem35,
 NOISEMATRIX_DATA.elem44 cscmatrix_elem44,
 NOISEMATRIX_DATA.elem45 cscmatrix_elem45,
 NOISEMATRIX_DATA.elem46 cscmatrix_elem46,
 NOISEMATRIX_DATA.elem55 cscmatrix_elem55,
 NOISEMATRIX_DATA.elem56 cscmatrix_elem56,
 NOISEMATRIX_DATA.elem57 cscmatrix_elem57,
 NOISEMATRIX_DATA.elem66 cscmatrix_elem66,
 NOISEMATRIX_DATA.elem67 cscmatrix_elem67,
 NOISEMATRIX_DATA.elem77 cscmatrix_elem77
FROM NOISEMATRIX_DATA,NOISEMATRIX_MAP
WHERE
 NOISEMATRIX_DATA.map_id=NOISEMATRIX_MAP.map_id
ORDER BY
 iov_value_id,
 map_id,
 vec_index;
