unit Unit1;

{$mode objfpc}{$H+}

interface

uses
  Classes, SysUtils, Forms, Controls, Graphics, Dialogs, StdCtrls, Math;

type

  { TForm1 }

  TForm1 = class(TForm)
    qEdit: TEdit;
    pEdit: TEdit;
    InputButton: TButton;
    ClearButton: TButton;
    Label1: TLabel;
    Label2: TLabel;
    HeaderListBox: TListBox;
    OpenDialog1: TOpenDialog;
    ORFPosButton: TButton;
    ORFNegButton: TButton;
    GroupBox1: TGroupBox;
    GroupBox2: TGroupBox;
    InputListBox: TListBox;
    OutputListBox: TListBox;
    procedure GroupBox1Click(Sender: TObject);
    procedure InputListBoxClick(Sender: TObject);
    procedure InputButtonClick(Sender: TObject);
    procedure ClearButtonClick(Sender: TObject);
    procedure Label2Click(Sender: TObject);
    procedure ORFPosButtonClick(Sender: TObject);
    procedure ORFNegButtonClick(Sender: TObject);
  private
    FSequence   : string;
    FTotalA, FTotalT, FTotalC, FTotalG : Integer;

    // HMM parameters (read from pEdit / qEdit)
    FP : Double;   // NC -> SC transition probability (start codon prior)
    FQ : Double;   // C3 -> ST transition probability (stop codon prior)

    // HMM log-space emission tables  [state 0..5][base 0..3]  (A=0 T=1 C=2 G=3)
    // States: 0=NC  1=SC  2=C1  3=C2  4=C3  5=ST
    FLogB : array[0..5, 0..3] of Double;

    // HMM log-space transition table  [from][to]
    FLogA : array[0..5, 0..5] of Double;

    // Log initial distribution
    FLogPi : array[0..5] of Double;

    procedure LoadSequence(const AFilename: string);
    procedure ComputeStats;
    procedure DisplayStats;
    function  GetP: Double;
    function  GetQ: Double;
    procedure InitHMM;
    function  ReverseComplement(const ASeq: string): string;
    procedure RunHMM(const AStrand: string; IsNegative: Boolean);
    procedure AddWrapped(const APrefix, ASeq: string);
  public

  end;

const
  NEG_INF = -1e300;

  // HMM state indices
  STATE_NC = 0;
  STATE_SC = 1;
  STATE_C1 = 2;
  STATE_C2 = 3;
  STATE_C3 = 4;
  STATE_ST = 5;

  MIN_ORF_LEN = 6;   // minimum nt beyond ATG before counting

var
  Form1: TForm1;

implementation

{$R *.lfm}

{ TForm1 }

// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────

function SafeLog(const X: Double): Double;
begin
  if X <= 0.0 then Result := NEG_INF
  else Result := Ln(X);
end;

function NtIndex(C: Char): Integer;
begin
  case UpCase(C) of
    'A': Result := 0;
    'T': Result := 1;
    'C': Result := 2;
    'G': Result := 3;
  else
    Result := -1;
  end;
end;

function TForm1.GetP: Double;
begin
  if not TryStrToFloat(Trim(pEdit.Text), Result) then
    Result := 0.001;
  if (Result <= 0) or (Result >= 1) then Result := 0.001;
end;

function TForm1.GetQ: Double;
begin
  if not TryStrToFloat(Trim(qEdit.Text), Result) then
    Result := 0.05;
  if (Result <= 0) or (Result >= 1) then Result := 0.05;
end;

// ─────────────────────────────────────────────────────────────
//  HMM initialisation
// ─────────────────────────────────────────────────────────────

procedure TForm1.InitHMM;
var
  I, J : Integer;
begin
  FP := GetP;
  FQ := GetQ;

  // Zero everything to NEG_INF
  for I := 0 to 5 do
  begin
    FLogPi[I] := NEG_INF;
    for J := 0 to 5 do FLogA[I][J] := NEG_INF;
    for J := 0 to 3 do FLogB[I][J] := NEG_INF;
  end;

  // ── Initial distribution: start in NC ──
  FLogPi[STATE_NC] := 0.0;   // ln(1.0)

  // ── Transitions ──
  FLogA[STATE_NC][STATE_NC] := SafeLog(1.0 - FP);
  FLogA[STATE_NC][STATE_SC] := SafeLog(FP);
  FLogA[STATE_SC][STATE_C1] := 0.0;   // deterministic
  FLogA[STATE_C1][STATE_C2] := 0.0;
  FLogA[STATE_C2][STATE_C3] := 0.0;
  FLogA[STATE_C3][STATE_C1] := SafeLog(1.0 - FQ);
  FLogA[STATE_C3][STATE_ST] := SafeLog(FQ);
  FLogA[STATE_ST][STATE_NC] := 0.0;   // deterministic

  // ── Emissions ──
  // NC: uniform background         A      T      C      G
  FLogB[STATE_NC][0] := SafeLog(0.25);
  FLogB[STATE_NC][1] := SafeLog(0.25);
  FLogB[STATE_NC][2] := SafeLog(0.25);
  FLogB[STATE_NC][3] := SafeLog(0.25);

  // SC: first base of ATG -> emits only A
  FLogB[STATE_SC][0] := 0.0;      // A  ln(1.0)
  FLogB[STATE_SC][1] := NEG_INF;  // T
  FLogB[STATE_SC][2] := NEG_INF;  // C
  FLogB[STATE_SC][3] := NEG_INF;  // G

  // C1: codon position 1 frequencies (typical phage/E.coli)
  FLogB[STATE_C1][0] := SafeLog(0.28);  // A
  FLogB[STATE_C1][1] := SafeLog(0.22);  // T
  FLogB[STATE_C1][2] := SafeLog(0.22);  // C
  FLogB[STATE_C1][3] := SafeLog(0.28);  // G

  // C2: codon position 2
  FLogB[STATE_C2][0] := SafeLog(0.28);
  FLogB[STATE_C2][1] := SafeLog(0.24);
  FLogB[STATE_C2][2] := SafeLog(0.20);
  FLogB[STATE_C2][3] := SafeLog(0.28);

  // C3: codon position 3
  FLogB[STATE_C3][0] := SafeLog(0.24);
  FLogB[STATE_C3][1] := SafeLog(0.28);
  FLogB[STATE_C3][2] := SafeLog(0.22);
  FLogB[STATE_C3][3] := SafeLog(0.26);

  // ST: averaged over TAA/TAG/TGA base composition
  FLogB[STATE_ST][0] := SafeLog(0.44);  // A
  FLogB[STATE_ST][1] := SafeLog(0.34);  // T
  FLogB[STATE_ST][2] := NEG_INF;        // C (never in stop codons)
  FLogB[STATE_ST][3] := SafeLog(0.22);  // G
end;

// ─────────────────────────────────────────────────────────────
//  FASTA loader
// ─────────────────────────────────────────────────────────────

procedure TForm1.LoadSequence(const AFilename: string);
var
  SL   : TStringList;
  I    : Integer;
  Line : string;
begin
  FSequence := '';
  InputListBox.Clear;
  HeaderListBox.Clear;

  SL := TStringList.Create;
  try
    SL.LoadFromFile(AFilename);
    for I := 0 to SL.Count - 1 do
    begin
      // Strip carriage returns for Windows-style files
      Line := Trim(StringReplace(SL[I], #13, '', [rfReplaceAll]));
      if Line = '' then Continue;

      // FASTA header
      if (Length(Line) > 0) and (Line[1] = '>') then
      begin
        HeaderListBox.Items.Add(Line);
        Continue;
      end;

      // Comment line
      if Line[1] = ';' then Continue;

      // Sequence data
      FSequence := FSequence + UpperCase(Line);
      InputListBox.Items.Add(UpperCase(Line));
    end;
  finally
    SL.Free;
  end;
end;

// ─────────────────────────────────────────────────────────────
//  Statistics
// ─────────────────────────────────────────────────────────────

procedure TForm1.ComputeStats;
var
  I  : Integer;
  Ch : Char;
begin
  FTotalA := 0; FTotalT := 0; FTotalC := 0; FTotalG := 0;
  for I := 1 to Length(FSequence) do
  begin
    Ch := FSequence[I];
    case Ch of
      'A': Inc(FTotalA);
      'T': Inc(FTotalT);
      'C': Inc(FTotalC);
      'G': Inc(FTotalG);
    end;
  end;
end;

procedure TForm1.DisplayStats;
var
  Total                       : Integer;
  FreqA, FreqT, FreqC, FreqG : Double;
begin
  Total := FTotalA + FTotalT + FTotalC + FTotalG;
  if Total > 0 then
  begin
    FreqA := FTotalA / Total;
    FreqT := FTotalT / Total;
    FreqC := FTotalC / Total;
    FreqG := FTotalG / Total;
  end
  else
  begin
    FreqA := 0; FreqT := 0; FreqC := 0; FreqG := 0;
  end;

  OutputListBox.Items.Add('Sequence Statistics');
  OutputListBox.Items.Add('Total Base Pairs : ' + IntToStr(Total));
  OutputListBox.Items.Add('total A : ' + IntToStr(FTotalA));
  OutputListBox.Items.Add('total T : ' + IntToStr(FTotalT));
  OutputListBox.Items.Add('total C : ' + IntToStr(FTotalC));
  OutputListBox.Items.Add('total G : ' + IntToStr(FTotalG));
  OutputListBox.Items.Add('freq A  : ' + FloatToStr(FreqA));
  OutputListBox.Items.Add('freq T  : ' + FloatToStr(FreqT));
  OutputListBox.Items.Add('freq C  : ' + FloatToStr(FreqC));
  OutputListBox.Items.Add('freq G  : ' + FloatToStr(FreqG));
  OutputListBox.Items.Add('accl freq : ' + FloatToStr(FreqA+FreqT+FreqC+FreqG));
  OutputListBox.Items.Add('');
  OutputListBox.Items.Add('HMM Parameters');
  OutputListBox.Items.Add('p (NC->SC, start prior) : ' + FloatToStr(FP));
  OutputListBox.Items.Add('q (C3->ST, stop  prior) : ' + FloatToStr(FQ));
  OutputListBox.Items.Add('Initial state           : NC (pi = 1.0)');
  OutputListBox.Items.Add('');
end;


function TForm1.ReverseComplement(const ASeq: string): string;
var
  I  : Integer;
  Ch : Char;
begin
  Result := '';
  for I := Length(ASeq) downto 1 do
  begin
    Ch := ASeq[I];
    case Ch of
      'A': Result := Result + 'T';
      'T': Result := Result + 'A';
      'C': Result := Result + 'G';
      'G': Result := Result + 'C';
    else
      Result := Result + Ch;
    end;
  end;
end;


procedure TForm1.AddWrapped(const APrefix, ASeq: string);
var
  WrapPos : Integer;
begin
  if Length(ASeq) <= 50 then
    OutputListBox.Items.Add(APrefix + ASeq)
  else
  begin
    OutputListBox.Items.Add(APrefix);
    WrapPos := 1;
    while WrapPos <= Length(ASeq) do
    begin
      OutputListBox.Items.Add('  ' + Copy(ASeq, WrapPos, 50));
      Inc(WrapPos, 50);
    end;
  end;
end;


function IsStopCodon(const ASeq: string; Pos: Integer): Boolean;
var
  Codon : string;
begin
  Codon := Copy(ASeq, Pos, 3);
  Result := (Codon = 'TAA') or (Codon = 'TAG') or (Codon = 'TGA');
end;


procedure TForm1.RunHMM(const AStrand: string; IsNegative: Boolean);
var
  SeqLen      : Integer;
  I, K        : Integer;
  ORFStart    : Integer;
  ORFSeq      : string;
  ORFCount    : Integer;
  FrameORFCount : array[0..2] of Integer;
  Frame       : Integer;
  FrameLabel  : string;
  StartPos, StopPos, ORFLen : Integer;

  CodingLogLik : Double;
  NCLogLik     : Double;
  NtIdx        : Integer;
  Codon        : string;
  CodonPos     : Integer;

  StateSeq     : array[0..2] of Integer;

begin
  SeqLen := Length(AStrand);
  ORFCount := 0;
  FrameORFCount[0] := 0;
  FrameORFCount[1] := 0;
  FrameORFCount[2] := 0;

  StateSeq[0] := STATE_C1;
  StateSeq[1] := STATE_C2;
  StateSeq[2] := STATE_C3;

  for Frame := 0 to 2 do
  begin
    if IsNegative then
      FrameLabel := '-' + IntToStr(Frame + 1)
    else
      if Frame = 0 then FrameLabel := '+0'
      else FrameLabel := '+' + IntToStr(Frame);

    I := 1 + Frame;

    while I + 2 <= SeqLen do
    begin
      if (AStrand[I] = 'A') and
         (AStrand[I+1] = 'T') and
         (AStrand[I+2] = 'G') then
      begin
        ORFStart := I;
        ORFSeq   := 'ATG';

        CodingLogLik := FLogB[STATE_SC][0]
                      + FLogB[STATE_C1][NtIndex('T')]
                      + FLogB[STATE_C2][NtIndex('G')];
        NCLogLik := FLogB[STATE_NC][0]
                  + FLogB[STATE_NC][1]
                  + FLogB[STATE_NC][3];
        K := I + 3;
        while K + 2 <= SeqLen do
        begin
          Codon := Copy(AStrand, K, 3);
          for CodonPos := 0 to 2 do
          begin
            NtIdx := NtIndex(Codon[CodonPos + 1]);
            if NtIdx >= 0 then
            begin
              CodingLogLik := CodingLogLik + FLogB[StateSeq[CodonPos]][NtIdx];
              NCLogLik     := NCLogLik     + FLogB[STATE_NC][NtIdx];
            end;
          end;

          ORFSeq := ORFSeq + Codon;
          if IsStopCodon(AStrand, K) then
          begin
            ORFLen := Length(ORFSeq);
            if (ORFLen >= MIN_ORF_LEN + 3) and (CodingLogLik > NCLogLik) then
            begin
              Inc(ORFCount);
              Inc(FrameORFCount[Frame]);
              if IsNegative then
              begin
                StartPos := SeqLen - (ORFStart + ORFLen - 1) + 1;
                StopPos  := SeqLen - ORFStart + 1;
              end
              else
              begin
                StartPos := ORFStart;
                StopPos  := ORFStart + ORFLen - 1;
              end;

              OutputListBox.Items.Add(
                Format('[Frame: %s][Start: %d  Stop: %d][Length: %d]' +
                       '[Score: %.2f]',
                  [FrameLabel, StartPos, StopPos, ORFLen,
                   CodingLogLik - NCLogLik])
              );
              AddWrapped('ORF[' + IntToStr(ORFCount) + ']: ', ORFSeq);
              OutputListBox.Items.Add('');
            end;

            K := K + 3;
            Break;
          end;

          K := K + 3;
        end;
       I := I + 3;
      end
      else
        Inc(I);
    end;
    OutputListBox.Items.Add(
      '[Frame ' + FrameLabel + '] HMM ORFs found: ' +
      IntToStr(FrameORFCount[Frame]) + ' ---'
    );
    OutputListBox.Items.Add('');
  end;
  OutputListBox.Items.Insert(
    OutputListBox.Items.IndexOf('') + 1,
    'HMM ORFs Found (total): ' + IntToStr(ORFCount)
  );
end;

procedure TForm1.InputButtonClick(Sender: TObject);
begin
  OpenDialog1.Filter :=
    'FASTA files (*.fasta;*.fa;*.txt)|*.fasta;*.fa;*.txt|All files (*.*)|*.*';
  OpenDialog1.Title := 'Open Sequence File';
  if OpenDialog1.Execute then
  begin
    LoadSequence(OpenDialog1.FileName);
    ComputeStats;
  end;
end;

procedure TForm1.ClearButtonClick(Sender: TObject);
begin
  OutputListBox.Clear;
  InputListBox.Clear;
  HeaderListBox.Clear;
  FSequence := '';
  FTotalA := 0; FTotalT := 0; FTotalC := 0; FTotalG := 0;
end;

procedure TForm1.Label2Click(Sender: TObject);
begin

end;

procedure TForm1.ORFPosButtonClick(Sender: TObject);
begin
  if FSequence = '' then
  begin
    ShowMessage('Please load a sequence file first.');
    Exit;
  end;
  OutputListBox.Clear;
  InitHMM;
  DisplayStats;
  RunHMM(FSequence, False);
end;

procedure TForm1.ORFNegButtonClick(Sender: TObject);
var
  RevComp : string;
begin
  if FSequence = '' then
  begin
    ShowMessage('Please load a sequence file first.');
    Exit;
  end;
  OutputListBox.Clear;
  InitHMM;
  DisplayStats;
  RevComp := ReverseComplement(FSequence);
  RunHMM(RevComp, True);
end;

procedure TForm1.InputListBoxClick(Sender: TObject);
begin

end;

procedure TForm1.GroupBox1Click(Sender: TObject);
begin

end;

end.
