# Step 1. 패키지 다운받기
- 줄리아의 패키지 형태는 docs/src/test/Readme.md로 구성된다.
- 패키지 소스를 특정 폴더에 다운받는다.
- 예) '/HDD/sugyeong/github/FP_test.jl' (FP_test.jl이라는 폴더이름에 저장함)

# Step 2. 패키지 설치
- (1) Julia 실행
- (2) ']' 눌러서 'pkg'실행
- (3) 다음 명령어를 입력한다.

```julia
dev /HDD/sugyeong/github/FP_test.jl
```
- result: 기존 'Juniper' 패키지가 다운받은 패키지로 바뀜

# Step 3. 패키지를 사용하는 .jl실행
- 수정된 패키지가 사용됨을 알 수 있음