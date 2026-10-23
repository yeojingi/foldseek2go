# Add MMseqs Command Skill

Use this skill when the user wants to add a new mmseqs command to this project (m2g).

## Files to modify (in order)

1. **Create** `src/util/<commandname>.cpp` (or `src/workflow/` for workflow commands, `src/func/` for function-related)
2. **Edit** `src/CommandDeclarations.h` — add forward declaration
3. **Edit** `src/MMseqsBase.cpp` — register in `baseCommands` vector
4. **Edit** `src/util/CMakeLists.txt` (or appropriate module) — add source file

커스텀 파라미터가 필요하면 추가로:

5. **Edit** `src/commons/Parameters.h` — 변수, `PARAMETER()` 매크로, 벡터 선언
6. **Edit** `src/commons/Parameters.cpp` — 생성자 초기화 목록 + 벡터 push_back

---

## Step 1: Implementation file

```cpp
#include "Parameters.h"
#include "Debug.h"

int <commandname>(int argc, const char **argv, const Command &command) {
    Parameters &par = Parameters::getInstance();
    par.parseParameters(argc, argv, command, true, 0, 0);

    // Implementation here
    // Input/output paths: par.db1, par.db2, par.db3, par.db4
    // Custom params: par.myParam

    return EXIT_SUCCESS;
}
```

## Step 2: Forward declaration

In `src/CommandDeclarations.h`, add:

```cpp
extern int <commandname>(int argc, const char **argv, const Command& command);
```

## Step 3: Register in baseCommands

In `src/MMseqsBase.cpp`, inside `std::vector<Command> baseCommands = { ... }`:

```cpp
{"<commandname>",  <commandname>,  &par.<paramgroup>,  COMMAND_DB,
        "Short description",
        "mmseqs <commandname> inputDB outputDB\n",
        "Author <email>",
        "<i:inputDB> <o:outputDB>",
        CITATION_MMSEQS2, {
            {"inputDB",  DbType::ACCESS_MODE_INPUT,  DbType::NEED_DATA, &DbValidator::sequenceDb},
            {"outputDB", DbType::ACCESS_MODE_OUTPUT, DbType::NEED_DATA, &DbValidator::sequenceDb}
        }},
```

## Step 4: Add to CMakeLists.txt

```cmake
util/<commandname>.cpp
```

---

## Adding a new parameter (커스텀 파라미터)

새 파라미터가 필요하면 **Parameters.h 2곳 + Parameters.cpp 2곳** 총 4곳을 수정합니다.

### Parameters.h — 수정 1: 값 저장 변수

클래스 멤버 변수 선언 영역에 추가 (관련 명령어 주석 근처에):

```cpp
// <commandname>
int    myParam = 0;        // 정수
float  myParam = 1.0f;     // 실수
bool   myParam = false;    // 플래그
std::string myParam = "";  // 문자열
```

### Parameters.h — 수정 2: PARAMETER 매크로 + 벡터 선언

`PARAMETER(...)` 선언 영역 (클래스 하단)에 추가:

```cpp
// <commandname>
PARAMETER(PARAM_MY_PARAM)
```

그리고 벡터 선언 영역에 전용 그룹 추가:

```cpp
std::vector<MMseqsParameter*> <commandname>;
```

### Parameters.cpp — 수정 3: 생성자 초기화 목록

`Parameters::Parameters():` 의 **초기화 목록 안**에 추가.
**반드시 여기에 넣어야 컴파일 에러가 나지 않음** — `MMseqsParameter`는 default constructor가 없기 때문에 `PARAMETER()` 로 선언한 모든 멤버는 반드시 초기화 목록에 있어야 함.

```cpp
// <commandname>
PARAM_MY_PARAM(PARAM_MY_PARAM_ID, "--my-param", "Display name",
    "Description of the option",
    typeid(int), (void *) &myParam,
    "^[0-1]{1}$",                        // 검증 regex
    MMseqsParameter::COMMAND_MISC),       // 표시 카테고리 (생략 가능)
```

기존 항목들과 쉼표(`,`)로 구분되어 있어야 함. 마지막 항목(`PARAM_HELP_LONG`)에는 쉼표 없음.

### Parameters.cpp — 수정 4: 벡터 push_back

생성자 **본문** (`{` 이후 블록)에서 관련 그룹 근처에 추가:

```cpp
// <commandname>
<commandname>.push_back(&PARAM_MY_PARAM);
<commandname>.push_back(&PARAM_THREADS);
<commandname>.push_back(&PARAM_V);
```

### MMseqsBase.cpp — 그룹 연결

`baseCommands`에서 해당 명령어의 `params` 필드를 전용 그룹으로 지정:

```cpp
{"<commandname>", <commandname>, &par.<commandname>, COMMAND_DB, ...}
//                                ^^^^^^^^^^^^^^^^^ par.<벡터이름>
```

---

## 파라미터 타입 & regex 참고

| 타입 | `typeid(...)` | regex 예시 |
|------|--------------|------------|
| 정수 | `typeid(int)` | `"^[0-9]{1}[0-9]*$"` |
| 0 또는 1 정수 | `typeid(int)` | `"^[0-1]{1}$"` |
| 실수 | `typeid(float)` | `"^[0-9]*(\\.[0-9]+)?$"` |
| 0.0-1.0 float | `typeid(float)` | `"^0(\\.[0-9]+)?|^1(\\.0+)?$"` |
| double | `typeid(double)` | `"^([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)"` |
| bool (flag) | `typeid(bool)` | `""` (regex 불필요) |
| 문자열 | `typeid(std::string)` | `""` |

## 파라미터 표시 카테고리 (`category` 필드)

| 상수 | 의미 |
|------|------|
| `COMMAND_COMMON` | 항상 표시 |
| `COMMAND_PREFILTER` | Prefilter 관련 |
| `COMMAND_ALIGN` | Alignment 관련 |
| `COMMAND_CLUST` | Clustering 관련 |
| `COMMAND_MISC` | 기타 (기본값, 생략 가능) |
| `COMMAND_EXPERT` | `--help` 기본 숨김 |
| `COMMAND_HIDDEN` | 항상 숨김 |

여러 조합: `MMseqsParameter::COMMAND_ALIGN | MMseqsParameter::COMMAND_EXPERT`

---

## Key reference

### Command categories (`mode` field)
| Flag | Category |
|------|----------|
| `COMMAND_MAIN` | Main workflows |
| `COMMAND_DB` | Database utilities |
| `COMMAND_FORMAT_CONVERSION` | Format conversion |
| `COMMAND_TAXONOMY` | Taxonomy |
| `COMMAND_ALIGNMENT` | Alignment |
| `COMMAND_CLUSTER` | Clustering |
| `COMMAND_HIDDEN` | Hidden from help |

### DB validators
| Validator | Accepts |
|-----------|---------|
| `&DbValidator::sequenceDb` | Sequence databases |
| `&DbValidator::alignmentDb` | Alignment results |
| `&DbValidator::clusterDb` | Cluster results |
| `&DbValidator::msaDb` | Multiple sequence alignments |
| `&DbValidator::flatfile` | FASTA / plain text |
| `&DbValidator::directory` | Directory path |

### 기존 parameter groups (전용 그룹 없이 재사용 가능)
| Group | 포함된 파라미터 |
|-------|----------------|
| `&par.onlyverbosity` | `-v` only |
| `&par.onlythreads` | `--threads`, `-v` |
| `&par.verbandcompression` | `--compressed`, `-v` |
| `&par.threadsandcompression` | `--threads`, `--compressed`, `-v` |
| `&par.align` | Alignment 관련 전체 |
| `&par.createdb` | DB 생성 관련 |

### Input/output access
| Constant | Meaning |
|----------|---------|
| `DbType::ACCESS_MODE_INPUT` | Read-only |
| `DbType::ACCESS_MODE_OUTPUT` | Write |
| `DbType::NEED_DATA` | Data file only |
| `DbType::NEED_HEADER` | Requires .h header file |
| `DbType::VARIADIC` | Variable number of inputs |

### Common parameter access
- `par.db1` ... `par.db4` — positional input/output paths (in `databases` vector order)
- `par.threads` — thread count
- `par.compressed` — output compression
- `par.sensitivity` — sensitivity setting

---

## Instructions for Claude

When this skill is invoked:
1. Ask the user: command name, source directory (`src/util/`, `src/workflow/`, `src/func/`), input/output DB types, brief description, and whether custom parameters are needed
2. Read `src/CommandDeclarations.h` to find the right insertion point
3. Read `src/MMseqsBase.cpp` to find a similar existing command as reference
4. **If custom parameters are needed:**
   - Read `src/commons/Parameters.h` — find the member variable section (near related command) and the `PARAMETER()` + vector declaration section (class bottom)
   - Read `src/commons/Parameters.cpp` — find the constructor initializer list (`:` to `{`) and the vector push_back section (constructor body)
   - Edit Parameters.h first (variables + PARAMETER macro + vector), then Parameters.cpp (initializer list entry + push_backs)
   - **Warn:** `PARAMETER()` 로 선언한 멤버는 반드시 초기화 목록에도 있어야 함. 빠지면 "no default constructor" 컴파일 에러 발생
5. Create the implementation file, then edit the remaining files in order
6. Remind the user to rebuild: `cd build && make`
