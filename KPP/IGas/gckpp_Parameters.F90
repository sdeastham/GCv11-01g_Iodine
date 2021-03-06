! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Parameter Module File
! 
! Generated by KPP-2.2.4_gc symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : gckpp_Parameters.f90
! Time                 : Sat Dec 24 18:09:14 2016
! Working directory    : /n/home13/seastham/GCStandard/Code/Code.v11-01g-Iodine/KPP/IGas
! Equation file        : gckpp.kpp
! Output root filename : gckpp
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE gckpp_Parameters

  USE gckpp_Precision
  PUBLIC
  SAVE


! NSPEC - Number of chemical species
  INTEGER, PARAMETER :: NSPEC = 809 
! NVAR - Number of Variable species
  INTEGER, PARAMETER :: NVAR = 801 
! NFLUX - Number of Reaction Flux species
  INTEGER, PARAMETER :: NFLUX = 616 
! NVARACT - Number of Active species
  INTEGER, PARAMETER :: NVARACT = 175 
! NFIX - Number of Fixed species
  INTEGER, PARAMETER :: NFIX = 8 
! NREACT - Number of reactions
  INTEGER, PARAMETER :: NREACT = 615 
! NVARST - Starting of variables in conc. vect.
  INTEGER, PARAMETER :: NVARST = 1 
! NFIXST - Starting of fixed in conc. vect.
  INTEGER, PARAMETER :: NFIXST = 802 
! NONZERO - Number of nonzero entries in Jacobian
  INTEGER, PARAMETER :: NONZERO = 3741 
! LU_NONZERO - Number of nonzero entries in LU factoriz. of Jacobian
  INTEGER, PARAMETER :: LU_NONZERO = 4302 
! CNVAR - (NVAR+1) Number of elements in compressed row format
  INTEGER, PARAMETER :: CNVAR = 802 
! NLOOKAT - Number of species to look at
  INTEGER, PARAMETER :: NLOOKAT = 0 
! NMONITOR - Number of species to monitor
  INTEGER, PARAMETER :: NMONITOR = 0 
! NMASS - Number of atoms to check mass balance
  INTEGER, PARAMETER :: NMASS = 1 

! Index declaration for variable species in C and VAR
!   VAR(ind_spc) = C(ind_spc)

  INTEGER, PARAMETER :: ind_MONX = 1 
  INTEGER, PARAMETER :: ind_CH2I2 = 2 
  INTEGER, PARAMETER :: ind_CH2ICl = 3 
  INTEGER, PARAMETER :: ind_CH2IBr = 4 
  INTEGER, PARAMETER :: ind_IBr = 5 
  INTEGER, PARAMETER :: ind_MSA = 6 
  INTEGER, PARAMETER :: ind_LISOPOH = 7 
  INTEGER, PARAMETER :: ind_LISOPNO3 = 8 
  INTEGER, PARAMETER :: ind_LBRO2H = 9 
  INTEGER, PARAMETER :: ind_LBRO2N = 10 
  INTEGER, PARAMETER :: ind_LTRO2H = 11 
  INTEGER, PARAMETER :: ind_LTRO2N = 12 
  INTEGER, PARAMETER :: ind_LXRO2H = 13 
  INTEGER, PARAMETER :: ind_LXRO2N = 14 
  INTEGER, PARAMETER :: ind_PYAC = 15 
  INTEGER, PARAMETER :: ind_CO2 = 16 
  INTEGER, PARAMETER :: ind_RR1 = 17 
  INTEGER, PARAMETER :: ind_RR2 = 18 
  INTEGER, PARAMETER :: ind_RR3 = 19 
  INTEGER, PARAMETER :: ind_RR4 = 20 
  INTEGER, PARAMETER :: ind_RR5 = 21 
  INTEGER, PARAMETER :: ind_RR6 = 22 
  INTEGER, PARAMETER :: ind_RR7 = 23 
  INTEGER, PARAMETER :: ind_RR8 = 24 
  INTEGER, PARAMETER :: ind_RR9 = 25 
  INTEGER, PARAMETER :: ind_RR10 = 26 
  INTEGER, PARAMETER :: ind_RR11 = 27 
  INTEGER, PARAMETER :: ind_RR12 = 28 
  INTEGER, PARAMETER :: ind_RR13 = 29 
  INTEGER, PARAMETER :: ind_RR14 = 30 
  INTEGER, PARAMETER :: ind_RR15 = 31 
  INTEGER, PARAMETER :: ind_RR16 = 32 
  INTEGER, PARAMETER :: ind_RR17 = 33 
  INTEGER, PARAMETER :: ind_RR18 = 34 
  INTEGER, PARAMETER :: ind_RR19 = 35 
  INTEGER, PARAMETER :: ind_RR20 = 36 
  INTEGER, PARAMETER :: ind_RR21 = 37 
  INTEGER, PARAMETER :: ind_RR22 = 38 
  INTEGER, PARAMETER :: ind_RR23 = 39 
  INTEGER, PARAMETER :: ind_RR24 = 40 
  INTEGER, PARAMETER :: ind_RR25 = 41 
  INTEGER, PARAMETER :: ind_RR26 = 42 
  INTEGER, PARAMETER :: ind_RR27 = 43 
  INTEGER, PARAMETER :: ind_RR28 = 44 
  INTEGER, PARAMETER :: ind_RR29 = 45 
  INTEGER, PARAMETER :: ind_RR30 = 46 
  INTEGER, PARAMETER :: ind_RR31 = 47 
  INTEGER, PARAMETER :: ind_RR32 = 48 
  INTEGER, PARAMETER :: ind_RR33 = 49 
  INTEGER, PARAMETER :: ind_RR34 = 50 
  INTEGER, PARAMETER :: ind_RR35 = 51 
  INTEGER, PARAMETER :: ind_RR36 = 52 
  INTEGER, PARAMETER :: ind_RR37 = 53 
  INTEGER, PARAMETER :: ind_RR38 = 54 
  INTEGER, PARAMETER :: ind_RR39 = 55 
  INTEGER, PARAMETER :: ind_RR40 = 56 
  INTEGER, PARAMETER :: ind_RR41 = 57 
  INTEGER, PARAMETER :: ind_RR42 = 58 
  INTEGER, PARAMETER :: ind_RR43 = 59 
  INTEGER, PARAMETER :: ind_RR44 = 60 
  INTEGER, PARAMETER :: ind_RR45 = 61 
  INTEGER, PARAMETER :: ind_RR46 = 62 
  INTEGER, PARAMETER :: ind_RR47 = 63 
  INTEGER, PARAMETER :: ind_RR48 = 64 
  INTEGER, PARAMETER :: ind_RR49 = 65 
  INTEGER, PARAMETER :: ind_RR50 = 66 
  INTEGER, PARAMETER :: ind_RR51 = 67 
  INTEGER, PARAMETER :: ind_RR52 = 68 
  INTEGER, PARAMETER :: ind_RR53 = 69 
  INTEGER, PARAMETER :: ind_RR54 = 70 
  INTEGER, PARAMETER :: ind_RR55 = 71 
  INTEGER, PARAMETER :: ind_RR56 = 72 
  INTEGER, PARAMETER :: ind_RR57 = 73 
  INTEGER, PARAMETER :: ind_RR58 = 74 
  INTEGER, PARAMETER :: ind_RR59 = 75 
  INTEGER, PARAMETER :: ind_RR60 = 76 
  INTEGER, PARAMETER :: ind_RR61 = 77 
  INTEGER, PARAMETER :: ind_RR62 = 78 
  INTEGER, PARAMETER :: ind_RR63 = 79 
  INTEGER, PARAMETER :: ind_RR64 = 80 
  INTEGER, PARAMETER :: ind_RR65 = 81 
  INTEGER, PARAMETER :: ind_RR66 = 82 
  INTEGER, PARAMETER :: ind_RR67 = 83 
  INTEGER, PARAMETER :: ind_RR68 = 84 
  INTEGER, PARAMETER :: ind_RR69 = 85 
  INTEGER, PARAMETER :: ind_RR70 = 86 
  INTEGER, PARAMETER :: ind_RR71 = 87 
  INTEGER, PARAMETER :: ind_RR72 = 88 
  INTEGER, PARAMETER :: ind_RR73 = 89 
  INTEGER, PARAMETER :: ind_RR74 = 90 
  INTEGER, PARAMETER :: ind_RR75 = 91 
  INTEGER, PARAMETER :: ind_RR76 = 92 
  INTEGER, PARAMETER :: ind_RR77 = 93 
  INTEGER, PARAMETER :: ind_RR78 = 94 
  INTEGER, PARAMETER :: ind_RR79 = 95 
  INTEGER, PARAMETER :: ind_RR80 = 96 
  INTEGER, PARAMETER :: ind_RR81 = 97 
  INTEGER, PARAMETER :: ind_RR82 = 98 
  INTEGER, PARAMETER :: ind_RR83 = 99 
  INTEGER, PARAMETER :: ind_RR84 = 100 
  INTEGER, PARAMETER :: ind_RR85 = 101 
  INTEGER, PARAMETER :: ind_RR86 = 102 
  INTEGER, PARAMETER :: ind_RR87 = 103 
  INTEGER, PARAMETER :: ind_RR88 = 104 
  INTEGER, PARAMETER :: ind_RR89 = 105 
  INTEGER, PARAMETER :: ind_RR90 = 106 
  INTEGER, PARAMETER :: ind_RR91 = 107 
  INTEGER, PARAMETER :: ind_RR92 = 108 
  INTEGER, PARAMETER :: ind_RR93 = 109 
  INTEGER, PARAMETER :: ind_RR94 = 110 
  INTEGER, PARAMETER :: ind_RR95 = 111 
  INTEGER, PARAMETER :: ind_RR96 = 112 
  INTEGER, PARAMETER :: ind_RR97 = 113 
  INTEGER, PARAMETER :: ind_RR98 = 114 
  INTEGER, PARAMETER :: ind_RR99 = 115 
  INTEGER, PARAMETER :: ind_RR100 = 116 
  INTEGER, PARAMETER :: ind_RR101 = 117 
  INTEGER, PARAMETER :: ind_RR102 = 118 
  INTEGER, PARAMETER :: ind_RR103 = 119 
  INTEGER, PARAMETER :: ind_RR104 = 120 
  INTEGER, PARAMETER :: ind_RR105 = 121 
  INTEGER, PARAMETER :: ind_RR106 = 122 
  INTEGER, PARAMETER :: ind_RR107 = 123 
  INTEGER, PARAMETER :: ind_RR108 = 124 
  INTEGER, PARAMETER :: ind_RR109 = 125 
  INTEGER, PARAMETER :: ind_RR110 = 126 
  INTEGER, PARAMETER :: ind_RR111 = 127 
  INTEGER, PARAMETER :: ind_RR112 = 128 
  INTEGER, PARAMETER :: ind_RR113 = 129 
  INTEGER, PARAMETER :: ind_RR114 = 130 
  INTEGER, PARAMETER :: ind_RR115 = 131 
  INTEGER, PARAMETER :: ind_RR116 = 132 
  INTEGER, PARAMETER :: ind_RR117 = 133 
  INTEGER, PARAMETER :: ind_RR118 = 134 
  INTEGER, PARAMETER :: ind_RR119 = 135 
  INTEGER, PARAMETER :: ind_RR120 = 136 
  INTEGER, PARAMETER :: ind_RR121 = 137 
  INTEGER, PARAMETER :: ind_RR122 = 138 
  INTEGER, PARAMETER :: ind_RR123 = 139 
  INTEGER, PARAMETER :: ind_RR124 = 140 
  INTEGER, PARAMETER :: ind_RR125 = 141 
  INTEGER, PARAMETER :: ind_RR126 = 142 
  INTEGER, PARAMETER :: ind_RR127 = 143 
  INTEGER, PARAMETER :: ind_RR128 = 144 
  INTEGER, PARAMETER :: ind_RR129 = 145 
  INTEGER, PARAMETER :: ind_RR130 = 146 
  INTEGER, PARAMETER :: ind_RR131 = 147 
  INTEGER, PARAMETER :: ind_RR132 = 148 
  INTEGER, PARAMETER :: ind_RR133 = 149 
  INTEGER, PARAMETER :: ind_RR134 = 150 
  INTEGER, PARAMETER :: ind_RR135 = 151 
  INTEGER, PARAMETER :: ind_RR136 = 152 
  INTEGER, PARAMETER :: ind_RR137 = 153 
  INTEGER, PARAMETER :: ind_RR138 = 154 
  INTEGER, PARAMETER :: ind_RR139 = 155 
  INTEGER, PARAMETER :: ind_RR140 = 156 
  INTEGER, PARAMETER :: ind_RR141 = 157 
  INTEGER, PARAMETER :: ind_RR142 = 158 
  INTEGER, PARAMETER :: ind_RR143 = 159 
  INTEGER, PARAMETER :: ind_RR144 = 160 
  INTEGER, PARAMETER :: ind_RR145 = 161 
  INTEGER, PARAMETER :: ind_RR146 = 162 
  INTEGER, PARAMETER :: ind_RR147 = 163 
  INTEGER, PARAMETER :: ind_RR148 = 164 
  INTEGER, PARAMETER :: ind_RR149 = 165 
  INTEGER, PARAMETER :: ind_RR150 = 166 
  INTEGER, PARAMETER :: ind_RR151 = 167 
  INTEGER, PARAMETER :: ind_RR152 = 168 
  INTEGER, PARAMETER :: ind_RR153 = 169 
  INTEGER, PARAMETER :: ind_RR154 = 170 
  INTEGER, PARAMETER :: ind_RR155 = 171 
  INTEGER, PARAMETER :: ind_RR156 = 172 
  INTEGER, PARAMETER :: ind_RR157 = 173 
  INTEGER, PARAMETER :: ind_RR158 = 174 
  INTEGER, PARAMETER :: ind_RR159 = 175 
  INTEGER, PARAMETER :: ind_RR160 = 176 
  INTEGER, PARAMETER :: ind_RR161 = 177 
  INTEGER, PARAMETER :: ind_RR162 = 178 
  INTEGER, PARAMETER :: ind_RR163 = 179 
  INTEGER, PARAMETER :: ind_RR164 = 180 
  INTEGER, PARAMETER :: ind_RR165 = 181 
  INTEGER, PARAMETER :: ind_RR166 = 182 
  INTEGER, PARAMETER :: ind_RR167 = 183 
  INTEGER, PARAMETER :: ind_RR168 = 184 
  INTEGER, PARAMETER :: ind_RR169 = 185 
  INTEGER, PARAMETER :: ind_RR170 = 186 
  INTEGER, PARAMETER :: ind_RR171 = 187 
  INTEGER, PARAMETER :: ind_RR172 = 188 
  INTEGER, PARAMETER :: ind_RR173 = 189 
  INTEGER, PARAMETER :: ind_RR174 = 190 
  INTEGER, PARAMETER :: ind_RR175 = 191 
  INTEGER, PARAMETER :: ind_RR176 = 192 
  INTEGER, PARAMETER :: ind_RR177 = 193 
  INTEGER, PARAMETER :: ind_RR178 = 194 
  INTEGER, PARAMETER :: ind_RR179 = 195 
  INTEGER, PARAMETER :: ind_RR180 = 196 
  INTEGER, PARAMETER :: ind_RR181 = 197 
  INTEGER, PARAMETER :: ind_RR182 = 198 
  INTEGER, PARAMETER :: ind_RR183 = 199 
  INTEGER, PARAMETER :: ind_RR184 = 200 
  INTEGER, PARAMETER :: ind_RR185 = 201 
  INTEGER, PARAMETER :: ind_RR186 = 202 
  INTEGER, PARAMETER :: ind_RR187 = 203 
  INTEGER, PARAMETER :: ind_RR188 = 204 
  INTEGER, PARAMETER :: ind_RR189 = 205 
  INTEGER, PARAMETER :: ind_RR190 = 206 
  INTEGER, PARAMETER :: ind_RR191 = 207 
  INTEGER, PARAMETER :: ind_RR192 = 208 
  INTEGER, PARAMETER :: ind_RR193 = 209 
  INTEGER, PARAMETER :: ind_RR194 = 210 
  INTEGER, PARAMETER :: ind_RR195 = 211 
  INTEGER, PARAMETER :: ind_RR196 = 212 
  INTEGER, PARAMETER :: ind_RR197 = 213 
  INTEGER, PARAMETER :: ind_RR198 = 214 
  INTEGER, PARAMETER :: ind_RR199 = 215 
  INTEGER, PARAMETER :: ind_RR200 = 216 
  INTEGER, PARAMETER :: ind_RR201 = 217 
  INTEGER, PARAMETER :: ind_RR202 = 218 
  INTEGER, PARAMETER :: ind_RR203 = 219 
  INTEGER, PARAMETER :: ind_RR204 = 220 
  INTEGER, PARAMETER :: ind_RR205 = 221 
  INTEGER, PARAMETER :: ind_RR206 = 222 
  INTEGER, PARAMETER :: ind_RR207 = 223 
  INTEGER, PARAMETER :: ind_RR208 = 224 
  INTEGER, PARAMETER :: ind_RR209 = 225 
  INTEGER, PARAMETER :: ind_RR210 = 226 
  INTEGER, PARAMETER :: ind_RR211 = 227 
  INTEGER, PARAMETER :: ind_RR212 = 228 
  INTEGER, PARAMETER :: ind_RR213 = 229 
  INTEGER, PARAMETER :: ind_RR214 = 230 
  INTEGER, PARAMETER :: ind_RR215 = 231 
  INTEGER, PARAMETER :: ind_RR216 = 232 
  INTEGER, PARAMETER :: ind_RR217 = 233 
  INTEGER, PARAMETER :: ind_RR218 = 234 
  INTEGER, PARAMETER :: ind_RR219 = 235 
  INTEGER, PARAMETER :: ind_RR220 = 236 
  INTEGER, PARAMETER :: ind_RR221 = 237 
  INTEGER, PARAMETER :: ind_RR222 = 238 
  INTEGER, PARAMETER :: ind_RR223 = 239 
  INTEGER, PARAMETER :: ind_RR224 = 240 
  INTEGER, PARAMETER :: ind_RR225 = 241 
  INTEGER, PARAMETER :: ind_RR226 = 242 
  INTEGER, PARAMETER :: ind_RR227 = 243 
  INTEGER, PARAMETER :: ind_RR228 = 244 
  INTEGER, PARAMETER :: ind_RR229 = 245 
  INTEGER, PARAMETER :: ind_RR230 = 246 
  INTEGER, PARAMETER :: ind_RR231 = 247 
  INTEGER, PARAMETER :: ind_RR232 = 248 
  INTEGER, PARAMETER :: ind_RR233 = 249 
  INTEGER, PARAMETER :: ind_RR234 = 250 
  INTEGER, PARAMETER :: ind_RR235 = 251 
  INTEGER, PARAMETER :: ind_RR236 = 252 
  INTEGER, PARAMETER :: ind_RR237 = 253 
  INTEGER, PARAMETER :: ind_RR238 = 254 
  INTEGER, PARAMETER :: ind_RR239 = 255 
  INTEGER, PARAMETER :: ind_RR240 = 256 
  INTEGER, PARAMETER :: ind_RR241 = 257 
  INTEGER, PARAMETER :: ind_RR242 = 258 
  INTEGER, PARAMETER :: ind_RR243 = 259 
  INTEGER, PARAMETER :: ind_RR244 = 260 
  INTEGER, PARAMETER :: ind_RR245 = 261 
  INTEGER, PARAMETER :: ind_RR246 = 262 
  INTEGER, PARAMETER :: ind_RR247 = 263 
  INTEGER, PARAMETER :: ind_RR248 = 264 
  INTEGER, PARAMETER :: ind_RR249 = 265 
  INTEGER, PARAMETER :: ind_RR250 = 266 
  INTEGER, PARAMETER :: ind_RR251 = 267 
  INTEGER, PARAMETER :: ind_RR252 = 268 
  INTEGER, PARAMETER :: ind_RR253 = 269 
  INTEGER, PARAMETER :: ind_RR254 = 270 
  INTEGER, PARAMETER :: ind_RR255 = 271 
  INTEGER, PARAMETER :: ind_RR256 = 272 
  INTEGER, PARAMETER :: ind_RR257 = 273 
  INTEGER, PARAMETER :: ind_RR258 = 274 
  INTEGER, PARAMETER :: ind_RR259 = 275 
  INTEGER, PARAMETER :: ind_RR260 = 276 
  INTEGER, PARAMETER :: ind_RR261 = 277 
  INTEGER, PARAMETER :: ind_RR262 = 278 
  INTEGER, PARAMETER :: ind_RR263 = 279 
  INTEGER, PARAMETER :: ind_RR264 = 280 
  INTEGER, PARAMETER :: ind_RR265 = 281 
  INTEGER, PARAMETER :: ind_RR266 = 282 
  INTEGER, PARAMETER :: ind_RR267 = 283 
  INTEGER, PARAMETER :: ind_RR268 = 284 
  INTEGER, PARAMETER :: ind_RR269 = 285 
  INTEGER, PARAMETER :: ind_RR270 = 286 
  INTEGER, PARAMETER :: ind_RR271 = 287 
  INTEGER, PARAMETER :: ind_RR272 = 288 
  INTEGER, PARAMETER :: ind_RR273 = 289 
  INTEGER, PARAMETER :: ind_RR274 = 290 
  INTEGER, PARAMETER :: ind_RR275 = 291 
  INTEGER, PARAMETER :: ind_RR276 = 292 
  INTEGER, PARAMETER :: ind_RR277 = 293 
  INTEGER, PARAMETER :: ind_RR278 = 294 
  INTEGER, PARAMETER :: ind_RR279 = 295 
  INTEGER, PARAMETER :: ind_RR280 = 296 
  INTEGER, PARAMETER :: ind_RR281 = 297 
  INTEGER, PARAMETER :: ind_RR282 = 298 
  INTEGER, PARAMETER :: ind_RR283 = 299 
  INTEGER, PARAMETER :: ind_RR284 = 300 
  INTEGER, PARAMETER :: ind_RR285 = 301 
  INTEGER, PARAMETER :: ind_RR286 = 302 
  INTEGER, PARAMETER :: ind_RR287 = 303 
  INTEGER, PARAMETER :: ind_RR288 = 304 
  INTEGER, PARAMETER :: ind_RR289 = 305 
  INTEGER, PARAMETER :: ind_RR290 = 306 
  INTEGER, PARAMETER :: ind_RR291 = 307 
  INTEGER, PARAMETER :: ind_RR292 = 308 
  INTEGER, PARAMETER :: ind_RR293 = 309 
  INTEGER, PARAMETER :: ind_RR294 = 310 
  INTEGER, PARAMETER :: ind_RR295 = 311 
  INTEGER, PARAMETER :: ind_RR296 = 312 
  INTEGER, PARAMETER :: ind_RR297 = 313 
  INTEGER, PARAMETER :: ind_RR298 = 314 
  INTEGER, PARAMETER :: ind_RR299 = 315 
  INTEGER, PARAMETER :: ind_RR300 = 316 
  INTEGER, PARAMETER :: ind_RR301 = 317 
  INTEGER, PARAMETER :: ind_RR302 = 318 
  INTEGER, PARAMETER :: ind_RR303 = 319 
  INTEGER, PARAMETER :: ind_RR304 = 320 
  INTEGER, PARAMETER :: ind_RR305 = 321 
  INTEGER, PARAMETER :: ind_RR306 = 322 
  INTEGER, PARAMETER :: ind_RR307 = 323 
  INTEGER, PARAMETER :: ind_RR308 = 324 
  INTEGER, PARAMETER :: ind_RR309 = 325 
  INTEGER, PARAMETER :: ind_RR310 = 326 
  INTEGER, PARAMETER :: ind_RR311 = 327 
  INTEGER, PARAMETER :: ind_RR312 = 328 
  INTEGER, PARAMETER :: ind_RR313 = 329 
  INTEGER, PARAMETER :: ind_RR314 = 330 
  INTEGER, PARAMETER :: ind_RR315 = 331 
  INTEGER, PARAMETER :: ind_RR316 = 332 
  INTEGER, PARAMETER :: ind_RR317 = 333 
  INTEGER, PARAMETER :: ind_RR318 = 334 
  INTEGER, PARAMETER :: ind_RR319 = 335 
  INTEGER, PARAMETER :: ind_RR320 = 336 
  INTEGER, PARAMETER :: ind_RR321 = 337 
  INTEGER, PARAMETER :: ind_RR322 = 338 
  INTEGER, PARAMETER :: ind_RR323 = 339 
  INTEGER, PARAMETER :: ind_RR324 = 340 
  INTEGER, PARAMETER :: ind_RR325 = 341 
  INTEGER, PARAMETER :: ind_RR326 = 342 
  INTEGER, PARAMETER :: ind_RR327 = 343 
  INTEGER, PARAMETER :: ind_RR328 = 344 
  INTEGER, PARAMETER :: ind_RR329 = 345 
  INTEGER, PARAMETER :: ind_RR330 = 346 
  INTEGER, PARAMETER :: ind_RR331 = 347 
  INTEGER, PARAMETER :: ind_RR332 = 348 
  INTEGER, PARAMETER :: ind_RR333 = 349 
  INTEGER, PARAMETER :: ind_RR334 = 350 
  INTEGER, PARAMETER :: ind_RR335 = 351 
  INTEGER, PARAMETER :: ind_RR336 = 352 
  INTEGER, PARAMETER :: ind_RR337 = 353 
  INTEGER, PARAMETER :: ind_RR338 = 354 
  INTEGER, PARAMETER :: ind_RR339 = 355 
  INTEGER, PARAMETER :: ind_RR340 = 356 
  INTEGER, PARAMETER :: ind_RR341 = 357 
  INTEGER, PARAMETER :: ind_RR342 = 358 
  INTEGER, PARAMETER :: ind_RR343 = 359 
  INTEGER, PARAMETER :: ind_RR344 = 360 
  INTEGER, PARAMETER :: ind_RR345 = 361 
  INTEGER, PARAMETER :: ind_RR346 = 362 
  INTEGER, PARAMETER :: ind_RR347 = 363 
  INTEGER, PARAMETER :: ind_RR348 = 364 
  INTEGER, PARAMETER :: ind_RR349 = 365 
  INTEGER, PARAMETER :: ind_RR350 = 366 
  INTEGER, PARAMETER :: ind_RR351 = 367 
  INTEGER, PARAMETER :: ind_RR352 = 368 
  INTEGER, PARAMETER :: ind_RR353 = 369 
  INTEGER, PARAMETER :: ind_RR354 = 370 
  INTEGER, PARAMETER :: ind_RR355 = 371 
  INTEGER, PARAMETER :: ind_RR356 = 372 
  INTEGER, PARAMETER :: ind_RR357 = 373 
  INTEGER, PARAMETER :: ind_RR358 = 374 
  INTEGER, PARAMETER :: ind_RR359 = 375 
  INTEGER, PARAMETER :: ind_RR360 = 376 
  INTEGER, PARAMETER :: ind_RR361 = 377 
  INTEGER, PARAMETER :: ind_RR362 = 378 
  INTEGER, PARAMETER :: ind_RR363 = 379 
  INTEGER, PARAMETER :: ind_RR364 = 380 
  INTEGER, PARAMETER :: ind_RR365 = 381 
  INTEGER, PARAMETER :: ind_RR366 = 382 
  INTEGER, PARAMETER :: ind_RR367 = 383 
  INTEGER, PARAMETER :: ind_RR368 = 384 
  INTEGER, PARAMETER :: ind_RR369 = 385 
  INTEGER, PARAMETER :: ind_RR370 = 386 
  INTEGER, PARAMETER :: ind_RR371 = 387 
  INTEGER, PARAMETER :: ind_RR372 = 388 
  INTEGER, PARAMETER :: ind_RR373 = 389 
  INTEGER, PARAMETER :: ind_RR374 = 390 
  INTEGER, PARAMETER :: ind_RR375 = 391 
  INTEGER, PARAMETER :: ind_RR376 = 392 
  INTEGER, PARAMETER :: ind_RR377 = 393 
  INTEGER, PARAMETER :: ind_RR378 = 394 
  INTEGER, PARAMETER :: ind_RR379 = 395 
  INTEGER, PARAMETER :: ind_RR380 = 396 
  INTEGER, PARAMETER :: ind_RR381 = 397 
  INTEGER, PARAMETER :: ind_RR382 = 398 
  INTEGER, PARAMETER :: ind_RR383 = 399 
  INTEGER, PARAMETER :: ind_RR384 = 400 
  INTEGER, PARAMETER :: ind_RR385 = 401 
  INTEGER, PARAMETER :: ind_RR386 = 402 
  INTEGER, PARAMETER :: ind_RR387 = 403 
  INTEGER, PARAMETER :: ind_RR388 = 404 
  INTEGER, PARAMETER :: ind_RR389 = 405 
  INTEGER, PARAMETER :: ind_RR390 = 406 
  INTEGER, PARAMETER :: ind_RR391 = 407 
  INTEGER, PARAMETER :: ind_RR392 = 408 
  INTEGER, PARAMETER :: ind_RR393 = 409 
  INTEGER, PARAMETER :: ind_RR394 = 410 
  INTEGER, PARAMETER :: ind_RR395 = 411 
  INTEGER, PARAMETER :: ind_RR396 = 412 
  INTEGER, PARAMETER :: ind_RR397 = 413 
  INTEGER, PARAMETER :: ind_RR398 = 414 
  INTEGER, PARAMETER :: ind_RR399 = 415 
  INTEGER, PARAMETER :: ind_RR400 = 416 
  INTEGER, PARAMETER :: ind_RR401 = 417 
  INTEGER, PARAMETER :: ind_RR402 = 418 
  INTEGER, PARAMETER :: ind_RR403 = 419 
  INTEGER, PARAMETER :: ind_RR404 = 420 
  INTEGER, PARAMETER :: ind_RR405 = 421 
  INTEGER, PARAMETER :: ind_RR406 = 422 
  INTEGER, PARAMETER :: ind_RR407 = 423 
  INTEGER, PARAMETER :: ind_RR408 = 424 
  INTEGER, PARAMETER :: ind_RR409 = 425 
  INTEGER, PARAMETER :: ind_RR410 = 426 
  INTEGER, PARAMETER :: ind_RR411 = 427 
  INTEGER, PARAMETER :: ind_RR412 = 428 
  INTEGER, PARAMETER :: ind_RR413 = 429 
  INTEGER, PARAMETER :: ind_RR414 = 430 
  INTEGER, PARAMETER :: ind_RR415 = 431 
  INTEGER, PARAMETER :: ind_RR416 = 432 
  INTEGER, PARAMETER :: ind_RR417 = 433 
  INTEGER, PARAMETER :: ind_RR418 = 434 
  INTEGER, PARAMETER :: ind_RR419 = 435 
  INTEGER, PARAMETER :: ind_RR420 = 436 
  INTEGER, PARAMETER :: ind_RR421 = 437 
  INTEGER, PARAMETER :: ind_RR422 = 438 
  INTEGER, PARAMETER :: ind_RR423 = 439 
  INTEGER, PARAMETER :: ind_RR424 = 440 
  INTEGER, PARAMETER :: ind_RR425 = 441 
  INTEGER, PARAMETER :: ind_RR426 = 442 
  INTEGER, PARAMETER :: ind_RR427 = 443 
  INTEGER, PARAMETER :: ind_RR428 = 444 
  INTEGER, PARAMETER :: ind_RR429 = 445 
  INTEGER, PARAMETER :: ind_RR430 = 446 
  INTEGER, PARAMETER :: ind_RR431 = 447 
  INTEGER, PARAMETER :: ind_RR432 = 448 
  INTEGER, PARAMETER :: ind_RR433 = 449 
  INTEGER, PARAMETER :: ind_RR434 = 450 
  INTEGER, PARAMETER :: ind_RR435 = 451 
  INTEGER, PARAMETER :: ind_RR436 = 452 
  INTEGER, PARAMETER :: ind_RR437 = 453 
  INTEGER, PARAMETER :: ind_RR438 = 454 
  INTEGER, PARAMETER :: ind_RR439 = 455 
  INTEGER, PARAMETER :: ind_RR440 = 456 
  INTEGER, PARAMETER :: ind_RR441 = 457 
  INTEGER, PARAMETER :: ind_RR442 = 458 
  INTEGER, PARAMETER :: ind_RR443 = 459 
  INTEGER, PARAMETER :: ind_RR444 = 460 
  INTEGER, PARAMETER :: ind_RR445 = 461 
  INTEGER, PARAMETER :: ind_RR446 = 462 
  INTEGER, PARAMETER :: ind_RR447 = 463 
  INTEGER, PARAMETER :: ind_RR448 = 464 
  INTEGER, PARAMETER :: ind_RR449 = 465 
  INTEGER, PARAMETER :: ind_RR450 = 466 
  INTEGER, PARAMETER :: ind_RR451 = 467 
  INTEGER, PARAMETER :: ind_RR452 = 468 
  INTEGER, PARAMETER :: ind_RR453 = 469 
  INTEGER, PARAMETER :: ind_RR454 = 470 
  INTEGER, PARAMETER :: ind_RR455 = 471 
  INTEGER, PARAMETER :: ind_RR456 = 472 
  INTEGER, PARAMETER :: ind_RR457 = 473 
  INTEGER, PARAMETER :: ind_RR458 = 474 
  INTEGER, PARAMETER :: ind_RR459 = 475 
  INTEGER, PARAMETER :: ind_RR460 = 476 
  INTEGER, PARAMETER :: ind_RR461 = 477 
  INTEGER, PARAMETER :: ind_RR462 = 478 
  INTEGER, PARAMETER :: ind_RR463 = 479 
  INTEGER, PARAMETER :: ind_RR464 = 480 
  INTEGER, PARAMETER :: ind_RR465 = 481 
  INTEGER, PARAMETER :: ind_RR466 = 482 
  INTEGER, PARAMETER :: ind_RR467 = 483 
  INTEGER, PARAMETER :: ind_RR468 = 484 
  INTEGER, PARAMETER :: ind_RR469 = 485 
  INTEGER, PARAMETER :: ind_RR470 = 486 
  INTEGER, PARAMETER :: ind_RR471 = 487 
  INTEGER, PARAMETER :: ind_RR472 = 488 
  INTEGER, PARAMETER :: ind_RR473 = 489 
  INTEGER, PARAMETER :: ind_RR474 = 490 
  INTEGER, PARAMETER :: ind_RR475 = 491 
  INTEGER, PARAMETER :: ind_RR476 = 492 
  INTEGER, PARAMETER :: ind_RR477 = 493 
  INTEGER, PARAMETER :: ind_RR478 = 494 
  INTEGER, PARAMETER :: ind_RR479 = 495 
  INTEGER, PARAMETER :: ind_RR480 = 496 
  INTEGER, PARAMETER :: ind_RR481 = 497 
  INTEGER, PARAMETER :: ind_RR482 = 498 
  INTEGER, PARAMETER :: ind_RR483 = 499 
  INTEGER, PARAMETER :: ind_RR484 = 500 
  INTEGER, PARAMETER :: ind_RR485 = 501 
  INTEGER, PARAMETER :: ind_RR486 = 502 
  INTEGER, PARAMETER :: ind_RR487 = 503 
  INTEGER, PARAMETER :: ind_RR488 = 504 
  INTEGER, PARAMETER :: ind_RR489 = 505 
  INTEGER, PARAMETER :: ind_RR490 = 506 
  INTEGER, PARAMETER :: ind_RR491 = 507 
  INTEGER, PARAMETER :: ind_RR492 = 508 
  INTEGER, PARAMETER :: ind_RR493 = 509 
  INTEGER, PARAMETER :: ind_RR494 = 510 
  INTEGER, PARAMETER :: ind_RR495 = 511 
  INTEGER, PARAMETER :: ind_RR496 = 512 
  INTEGER, PARAMETER :: ind_RR497 = 513 
  INTEGER, PARAMETER :: ind_RR498 = 514 
  INTEGER, PARAMETER :: ind_RR499 = 515 
  INTEGER, PARAMETER :: ind_RR500 = 516 
  INTEGER, PARAMETER :: ind_RR501 = 517 
  INTEGER, PARAMETER :: ind_RR502 = 518 
  INTEGER, PARAMETER :: ind_RR503 = 519 
  INTEGER, PARAMETER :: ind_RR504 = 520 
  INTEGER, PARAMETER :: ind_RR505 = 521 
  INTEGER, PARAMETER :: ind_RR506 = 522 
  INTEGER, PARAMETER :: ind_RR507 = 523 
  INTEGER, PARAMETER :: ind_RR508 = 524 
  INTEGER, PARAMETER :: ind_RR509 = 525 
  INTEGER, PARAMETER :: ind_RR510 = 526 
  INTEGER, PARAMETER :: ind_RR511 = 527 
  INTEGER, PARAMETER :: ind_RR512 = 528 
  INTEGER, PARAMETER :: ind_RR513 = 529 
  INTEGER, PARAMETER :: ind_RR514 = 530 
  INTEGER, PARAMETER :: ind_RR515 = 531 
  INTEGER, PARAMETER :: ind_RR516 = 532 
  INTEGER, PARAMETER :: ind_RR517 = 533 
  INTEGER, PARAMETER :: ind_RR518 = 534 
  INTEGER, PARAMETER :: ind_RR519 = 535 
  INTEGER, PARAMETER :: ind_RR520 = 536 
  INTEGER, PARAMETER :: ind_RR521 = 537 
  INTEGER, PARAMETER :: ind_RR522 = 538 
  INTEGER, PARAMETER :: ind_RR523 = 539 
  INTEGER, PARAMETER :: ind_RR524 = 540 
  INTEGER, PARAMETER :: ind_RR525 = 541 
  INTEGER, PARAMETER :: ind_RR526 = 542 
  INTEGER, PARAMETER :: ind_RR527 = 543 
  INTEGER, PARAMETER :: ind_RR528 = 544 
  INTEGER, PARAMETER :: ind_RR529 = 545 
  INTEGER, PARAMETER :: ind_RR530 = 546 
  INTEGER, PARAMETER :: ind_RR531 = 547 
  INTEGER, PARAMETER :: ind_RR532 = 548 
  INTEGER, PARAMETER :: ind_RR533 = 549 
  INTEGER, PARAMETER :: ind_RR534 = 550 
  INTEGER, PARAMETER :: ind_RR535 = 551 
  INTEGER, PARAMETER :: ind_RR536 = 552 
  INTEGER, PARAMETER :: ind_RR537 = 553 
  INTEGER, PARAMETER :: ind_RR538 = 554 
  INTEGER, PARAMETER :: ind_RR539 = 555 
  INTEGER, PARAMETER :: ind_RR540 = 556 
  INTEGER, PARAMETER :: ind_RR541 = 557 
  INTEGER, PARAMETER :: ind_RR542 = 558 
  INTEGER, PARAMETER :: ind_RR543 = 559 
  INTEGER, PARAMETER :: ind_RR544 = 560 
  INTEGER, PARAMETER :: ind_RR545 = 561 
  INTEGER, PARAMETER :: ind_RR546 = 562 
  INTEGER, PARAMETER :: ind_RR547 = 563 
  INTEGER, PARAMETER :: ind_RR548 = 564 
  INTEGER, PARAMETER :: ind_RR549 = 565 
  INTEGER, PARAMETER :: ind_RR550 = 566 
  INTEGER, PARAMETER :: ind_RR551 = 567 
  INTEGER, PARAMETER :: ind_RR552 = 568 
  INTEGER, PARAMETER :: ind_RR553 = 569 
  INTEGER, PARAMETER :: ind_RR554 = 570 
  INTEGER, PARAMETER :: ind_RR555 = 571 
  INTEGER, PARAMETER :: ind_RR556 = 572 
  INTEGER, PARAMETER :: ind_RR557 = 573 
  INTEGER, PARAMETER :: ind_RR558 = 574 
  INTEGER, PARAMETER :: ind_RR559 = 575 
  INTEGER, PARAMETER :: ind_RR560 = 576 
  INTEGER, PARAMETER :: ind_RR561 = 577 
  INTEGER, PARAMETER :: ind_RR562 = 578 
  INTEGER, PARAMETER :: ind_RR563 = 579 
  INTEGER, PARAMETER :: ind_RR564 = 580 
  INTEGER, PARAMETER :: ind_RR565 = 581 
  INTEGER, PARAMETER :: ind_RR566 = 582 
  INTEGER, PARAMETER :: ind_RR567 = 583 
  INTEGER, PARAMETER :: ind_RR568 = 584 
  INTEGER, PARAMETER :: ind_RR569 = 585 
  INTEGER, PARAMETER :: ind_RR570 = 586 
  INTEGER, PARAMETER :: ind_RR571 = 587 
  INTEGER, PARAMETER :: ind_RR572 = 588 
  INTEGER, PARAMETER :: ind_RR573 = 589 
  INTEGER, PARAMETER :: ind_RR574 = 590 
  INTEGER, PARAMETER :: ind_RR575 = 591 
  INTEGER, PARAMETER :: ind_RR576 = 592 
  INTEGER, PARAMETER :: ind_RR577 = 593 
  INTEGER, PARAMETER :: ind_RR578 = 594 
  INTEGER, PARAMETER :: ind_RR579 = 595 
  INTEGER, PARAMETER :: ind_RR580 = 596 
  INTEGER, PARAMETER :: ind_RR581 = 597 
  INTEGER, PARAMETER :: ind_RR582 = 598 
  INTEGER, PARAMETER :: ind_RR583 = 599 
  INTEGER, PARAMETER :: ind_RR584 = 600 
  INTEGER, PARAMETER :: ind_RR585 = 601 
  INTEGER, PARAMETER :: ind_RR586 = 602 
  INTEGER, PARAMETER :: ind_RR587 = 603 
  INTEGER, PARAMETER :: ind_RR588 = 604 
  INTEGER, PARAMETER :: ind_RR589 = 605 
  INTEGER, PARAMETER :: ind_RR590 = 606 
  INTEGER, PARAMETER :: ind_RR591 = 607 
  INTEGER, PARAMETER :: ind_RR592 = 608 
  INTEGER, PARAMETER :: ind_RR593 = 609 
  INTEGER, PARAMETER :: ind_RR594 = 610 
  INTEGER, PARAMETER :: ind_RR595 = 611 
  INTEGER, PARAMETER :: ind_RR596 = 612 
  INTEGER, PARAMETER :: ind_RR597 = 613 
  INTEGER, PARAMETER :: ind_RR598 = 614 
  INTEGER, PARAMETER :: ind_RR599 = 615 
  INTEGER, PARAMETER :: ind_RR600 = 616 
  INTEGER, PARAMETER :: ind_RR601 = 617 
  INTEGER, PARAMETER :: ind_RR602 = 618 
  INTEGER, PARAMETER :: ind_RR603 = 619 
  INTEGER, PARAMETER :: ind_RR604 = 620 
  INTEGER, PARAMETER :: ind_RR605 = 621 
  INTEGER, PARAMETER :: ind_RR606 = 622 
  INTEGER, PARAMETER :: ind_RR607 = 623 
  INTEGER, PARAMETER :: ind_RR608 = 624 
  INTEGER, PARAMETER :: ind_RR609 = 625 
  INTEGER, PARAMETER :: ind_RR610 = 626 
  INTEGER, PARAMETER :: ind_RR611 = 627 
  INTEGER, PARAMETER :: ind_RR612 = 628 
  INTEGER, PARAMETER :: ind_RR613 = 629 
  INTEGER, PARAMETER :: ind_RR614 = 630 
  INTEGER, PARAMETER :: ind_RR615 = 631 
  INTEGER, PARAMETER :: ind_I2O4 = 632 
  INTEGER, PARAMETER :: ind_BENZ = 633 
  INTEGER, PARAMETER :: ind_TOLU = 634 
  INTEGER, PARAMETER :: ind_XYLE = 635 
  INTEGER, PARAMETER :: ind_CH3CCl3 = 636 
  INTEGER, PARAMETER :: ind_I2O2 = 637 
  INTEGER, PARAMETER :: ind_CCl4 = 638 
  INTEGER, PARAMETER :: ind_CFC11 = 639 
  INTEGER, PARAMETER :: ind_CFC12 = 640 
  INTEGER, PARAMETER :: ind_CFC113 = 641 
  INTEGER, PARAMETER :: ind_CFC114 = 642 
  INTEGER, PARAMETER :: ind_CFC115 = 643 
  INTEGER, PARAMETER :: ind_H1301 = 644 
  INTEGER, PARAMETER :: ind_H2402 = 645 
  INTEGER, PARAMETER :: ind_CH3I = 646 
  INTEGER, PARAMETER :: ind_ICl = 647 
  INTEGER, PARAMETER :: ind_I2O3 = 648 
  INTEGER, PARAMETER :: ind_PPN = 649 
  INTEGER, PARAMETER :: ind_BrNO2 = 650 
  INTEGER, PARAMETER :: ind_IEPOX = 651 
  INTEGER, PARAMETER :: ind_PMNN = 652 
  INTEGER, PARAMETER :: ind_H1211 = 653 
  INTEGER, PARAMETER :: ind_INO = 654 
  INTEGER, PARAMETER :: ind_IONO = 655 
  INTEGER, PARAMETER :: ind_N2O = 656 
  INTEGER, PARAMETER :: ind_BRO2 = 657 
  INTEGER, PARAMETER :: ind_TRO2 = 658 
  INTEGER, PARAMETER :: ind_XRO2 = 659 
  INTEGER, PARAMETER :: ind_OCS = 660 
  INTEGER, PARAMETER :: ind_N = 661 
  INTEGER, PARAMETER :: ind_PAN = 662 
  INTEGER, PARAMETER :: ind_HI = 663 
  INTEGER, PARAMETER :: ind_MAP = 664 
  INTEGER, PARAMETER :: ind_Cl2O2 = 665 
  INTEGER, PARAMETER :: ind_CHCl3 = 666 
  INTEGER, PARAMETER :: ind_CH2Cl2 = 667 
  INTEGER, PARAMETER :: ind_CHBr3 = 668 
  INTEGER, PARAMETER :: ind_CH2Br2 = 669 
  INTEGER, PARAMETER :: ind_MPN = 670 
  INTEGER, PARAMETER :: ind_OIO = 671 
  INTEGER, PARAMETER :: ind_HCFC123 = 672 
  INTEGER, PARAMETER :: ind_ETP = 673 
  INTEGER, PARAMETER :: ind_HNO2 = 674 
  INTEGER, PARAMETER :: ind_ClNO2 = 675 
  INTEGER, PARAMETER :: ind_HCFC141b = 676 
  INTEGER, PARAMETER :: ind_HCFC142b = 677 
  INTEGER, PARAMETER :: ind_RA3P = 678 
  INTEGER, PARAMETER :: ind_RB3P = 679 
  INTEGER, PARAMETER :: ind_DMS = 680 
  INTEGER, PARAMETER :: ind_HCFC22 = 681 
  INTEGER, PARAMETER :: ind_CH3Br = 682 
  INTEGER, PARAMETER :: ind_CH3Cl = 683 
  INTEGER, PARAMETER :: ind_HNO4 = 684 
  INTEGER, PARAMETER :: ind_MAOP = 685 
  INTEGER, PARAMETER :: ind_ClOO = 686 
  INTEGER, PARAMETER :: ind_HOI = 687 
  INTEGER, PARAMETER :: ind_RP = 688 
  INTEGER, PARAMETER :: ind_OClO = 689 
  INTEGER, PARAMETER :: ind_PP = 690 
  INTEGER, PARAMETER :: ind_PRPN = 691 
  INTEGER, PARAMETER :: ind_SO4 = 692 
  INTEGER, PARAMETER :: ind_ETHLN = 693 
  INTEGER, PARAMETER :: ind_ALK4 = 694 
  INTEGER, PARAMETER :: ind_R4P = 695 
  INTEGER, PARAMETER :: ind_BrCl = 696 
  INTEGER, PARAMETER :: ind_IAP = 697 
  INTEGER, PARAMETER :: ind_RIP = 698 
  INTEGER, PARAMETER :: ind_VRP = 699 
  INTEGER, PARAMETER :: ind_C3H8 = 700 
  INTEGER, PARAMETER :: ind_Br2 = 701 
  INTEGER, PARAMETER :: ind_ATOOH = 702 
  INTEGER, PARAMETER :: ind_IONO2 = 703 
  INTEGER, PARAMETER :: ind_DHMOB = 704 
  INTEGER, PARAMETER :: ind_MOBA = 705 
  INTEGER, PARAMETER :: ind_MP = 706 
  INTEGER, PARAMETER :: ind_BrSALA = 707 
  INTEGER, PARAMETER :: ind_BrSALC = 708 
  INTEGER, PARAMETER :: ind_MRP = 709 
  INTEGER, PARAMETER :: ind_N2O5 = 710 
  INTEGER, PARAMETER :: ind_ISNOHOO = 711 
  INTEGER, PARAMETER :: ind_ISNP = 712 
  INTEGER, PARAMETER :: ind_I2 = 713 
  INTEGER, PARAMETER :: ind_ISOPNB = 714 
  INTEGER, PARAMETER :: ind_C2H6 = 715 
  INTEGER, PARAMETER :: ind_IEPOXOO = 716 
  INTEGER, PARAMETER :: ind_MACRNO2 = 717 
  INTEGER, PARAMETER :: ind_ROH = 718 
  INTEGER, PARAMETER :: ind_MOBAOO = 719 
  INTEGER, PARAMETER :: ind_DIBOO = 720 
  INTEGER, PARAMETER :: ind_ISNOOB = 721 
  INTEGER, PARAMETER :: ind_PMN = 722 
  INTEGER, PARAMETER :: ind_INPN = 723 
  INTEGER, PARAMETER :: ind_MVKN = 724 
  INTEGER, PARAMETER :: ind_BrNO3 = 725 
  INTEGER, PARAMETER :: ind_H = 726 
  INTEGER, PARAMETER :: ind_Cl2 = 727 
  INTEGER, PARAMETER :: ind_CH4 = 728 
  INTEGER, PARAMETER :: ind_ISOPND = 729 
  INTEGER, PARAMETER :: ind_MVKOO = 730 
  INTEGER, PARAMETER :: ind_CH3CHOO = 731 
  INTEGER, PARAMETER :: ind_GLYX = 732 
  INTEGER, PARAMETER :: ind_MAOPO2 = 733 
  INTEGER, PARAMETER :: ind_PROPNN = 734 
  INTEGER, PARAMETER :: ind_GAOO = 735 
  INTEGER, PARAMETER :: ind_MGLYOO = 736 
  INTEGER, PARAMETER :: ind_A3O2 = 737 
  INTEGER, PARAMETER :: ind_MACRN = 738 
  INTEGER, PARAMETER :: ind_CH2OO = 739 
  INTEGER, PARAMETER :: ind_PRN1 = 740 
  INTEGER, PARAMETER :: ind_MGLOO = 741 
  INTEGER, PARAMETER :: ind_PO2 = 742 
  INTEGER, PARAMETER :: ind_B3O2 = 743 
  INTEGER, PARAMETER :: ind_ISNOOA = 744 
  INTEGER, PARAMETER :: ind_MAN2 = 745 
  INTEGER, PARAMETER :: ind_ISOP = 746 
  INTEGER, PARAMETER :: ind_ACET = 747 
  INTEGER, PARAMETER :: ind_H2O2 = 748 
  INTEGER, PARAMETER :: ind_KO2 = 749 
  INTEGER, PARAMETER :: ind_I = 750 
  INTEGER, PARAMETER :: ind_GLYC = 751 
  INTEGER, PARAMETER :: ind_HC5OO = 752 
  INTEGER, PARAMETER :: ind_MRO2 = 753 
  INTEGER, PARAMETER :: ind_RCO3 = 754 
  INTEGER, PARAMETER :: ind_R4N1 = 755 
  INTEGER, PARAMETER :: ind_VRO2 = 756 
  INTEGER, PARAMETER :: ind_IO = 757 
  INTEGER, PARAMETER :: ind_HOBr = 758 
  INTEGER, PARAMETER :: ind_HNO3 = 759 
  INTEGER, PARAMETER :: ind_PRPE = 760 
  INTEGER, PARAMETER :: ind_ATO2 = 761 
  INTEGER, PARAMETER :: ind_ISN1 = 762 
  INTEGER, PARAMETER :: ind_MAO3 = 763 
  INTEGER, PARAMETER :: ind_HOCl = 764 
  INTEGER, PARAMETER :: ind_ClNO3 = 765 
  INTEGER, PARAMETER :: ind_ETO2 = 766 
  INTEGER, PARAMETER :: ind_HAC = 767 
  INTEGER, PARAMETER :: ind_RIO2 = 768 
  INTEGER, PARAMETER :: ind_ISOPNBO2 = 769 
  INTEGER, PARAMETER :: ind_MGLY = 770 
  INTEGER, PARAMETER :: ind_ISOPNDO2 = 771 
  INTEGER, PARAMETER :: ind_MACROO = 772 
  INTEGER, PARAMETER :: ind_HC5 = 773 
  INTEGER, PARAMETER :: ind_R4O2 = 774 
  INTEGER, PARAMETER :: ind_INO2 = 775 
  INTEGER, PARAMETER :: ind_R4N2 = 776 
  INTEGER, PARAMETER :: ind_RCHO = 777 
  INTEGER, PARAMETER :: ind_BrO = 778 
  INTEGER, PARAMETER :: ind_CH2O = 779 
  INTEGER, PARAMETER :: ind_MACR = 780 
  INTEGER, PARAMETER :: ind_ALD2 = 781 
  INTEGER, PARAMETER :: ind_MVK = 782 
  INTEGER, PARAMETER :: ind_MEK = 783 
  INTEGER, PARAMETER :: ind_MCO3 = 784 
  INTEGER, PARAMETER :: ind_SO2 = 785 
  INTEGER, PARAMETER :: ind_HBr = 786 
  INTEGER, PARAMETER :: ind_HCl = 787 
  INTEGER, PARAMETER :: ind_NO = 788 
  INTEGER, PARAMETER :: ind_O1D = 789 
  INTEGER, PARAMETER :: ind_Br = 790 
  INTEGER, PARAMETER :: ind_NO2 = 791 
  INTEGER, PARAMETER :: ind_CO = 792 
  INTEGER, PARAMETER :: ind_O3 = 793 
  INTEGER, PARAMETER :: ind_MO2 = 794 
  INTEGER, PARAMETER :: ind_NO3 = 795 
  INTEGER, PARAMETER :: ind_H2O = 796 
  INTEGER, PARAMETER :: ind_O = 797 
  INTEGER, PARAMETER :: ind_ClO = 798 
  INTEGER, PARAMETER :: ind_HO2 = 799 
  INTEGER, PARAMETER :: ind_OH = 800 
  INTEGER, PARAMETER :: ind_Cl = 801 

! Index declaration for fixed species in C
!   C(ind_spc)

  INTEGER, PARAMETER :: ind_ACTA = 802 
  INTEGER, PARAMETER :: ind_EOH = 803 
  INTEGER, PARAMETER :: ind_H2 = 804 
  INTEGER, PARAMETER :: ind_HCOOH = 805 
  INTEGER, PARAMETER :: ind_MOH = 806 
  INTEGER, PARAMETER :: ind_N2 = 807 
  INTEGER, PARAMETER :: ind_O2 = 808 
  INTEGER, PARAMETER :: ind_RCOOH = 809 

! Index declaration for fixed species in FIX
!    FIX(indf_spc) = C(ind_spc) = C(NVAR+indf_spc)

  INTEGER, PARAMETER :: indf_ACTA = 1 
  INTEGER, PARAMETER :: indf_EOH = 2 
  INTEGER, PARAMETER :: indf_H2 = 3 
  INTEGER, PARAMETER :: indf_HCOOH = 4 
  INTEGER, PARAMETER :: indf_MOH = 5 
  INTEGER, PARAMETER :: indf_N2 = 6 
  INTEGER, PARAMETER :: indf_O2 = 7 
  INTEGER, PARAMETER :: indf_RCOOH = 8 

END MODULE gckpp_Parameters

